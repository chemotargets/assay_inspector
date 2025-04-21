#! /usr/bin/env python

__author__ = "Luca Menestrina"
__date__ = "20250130"
__copyright__ = "Copyright 2025, Chemotargets"
__license__ = ""
__credits__ = ["Data Science & Translational Research Group", "Ismael Gomez", "Ricard Garcia"]
__maintainer__ = "Luca Menestrina"
__version__ = "20250224"
__deprecated__ = False

import os
from rdkit import Chem
import pandas as pd
import numpy as np
from tqdm import tqdm
from copy import deepcopy

from AI_Utils import parallel_apply, ray_manager, logging, molFromSmiles, standardize
from AI_MoleculeInfo import MoleculeInfo

N_CPUS = 24  # TODO: move to configs
TASK_REGRESSION = "REGRESSION"  # TODO: move to configs
TASK_CLASSIFICATION = "CLASSIFICATION"  # TODO: move to configs
VALUES_ROUNDING = 6  # TODO: move to configs

###     MoleculeData class

class MoleculeData():

    """
    Class to handle groups of molecules as a whole (including transforming them into the formats required by DataForModeling).

    This class can be represented as a dataframe, a dictionary of MoleculeInfo objects, or both. Data can be loaded from files or from database,
    depending on the needs of the application.
    """
    ###      Define constants

    ##     Molecule data-related names
    NAME_INCHIKEY = "inchikey"
    NAME_SMILES_RAW = "smiles"
    NAME_SMILES_STD = "smiles_std"
    NAME_MOLOBJ = "molObj"
    NAME_VALUE = "value"
    NAME_ENDPOINT = "endpoint"
    NAME_REF = "ref"

    ##      Item names  # TODO: check names
    NAME_ID = "id"
    NAME_ITEM = NAME_SMILES_RAW  # TODO: it definitely needs a better name

    ##      Mandatory columns
    MANDATORY_COLS = [NAME_ID, NAME_SMILES_RAW, NAME_SMILES_STD, NAME_INCHIKEY]

    ###
    def __init__(self,
                 source = None,
                 ids = None,
                 standardized = False,
                 deduplicated = False,
                 preprocessed = False,
                 ):
        
        self._ids = ids
        self.is_standardized = standardized
        self.is_deduplicated = deduplicated
        self.is_preprocessed = preprocessed
        self.featuresInDF_list = []
        self._dataframe = None

        ##      Store file path
        self.source = source

        if not source is None:  # TODO: separate into loadData (also standardizes the DF)
            if isinstance(self.source, str):
                self.__loadFile()
            elif isinstance(self.source, pd.DataFrame):
                self._dataframe = deepcopy(self.source)
            else:
                logging.critical("Unknown source format")
            
            ##      Associate ids to items
            if self._ids is None:
                if self.NAME_ID in self._dataframe.columns:
                    self._ids = self._dataframe[self.NAME_ID].to_list()
                else:
                    self._ids = range(len(self._dataframe))
            if self.NAME_ID not in self._dataframe.columns:
                self._dataframe[self.NAME_ID] = self._ids

            ##      Standardization of data
            if not self.is_standardized:
                self.__standardize()

            ##      Remove not finite endpoint values
            self._dataframe = self._dataframe.replace((np.inf, -np.inf), np.nan)
            if self.NAME_VALUE in self._dataframe.columns:
                self._dataframe = self._dataframe.dropna(subset=[self.NAME_VALUE])

        ##      Initialize dictionary of MoleculeInfo objects
        self._molecules = None

    def __loadFile(self):
        if self.source is None:
            logging.warning(f"Trying to load file with no path")
            return
        
        ##      Check that file exists
        if not os.path.isfile(self.source):
            logging.error(f"File {self.source} doesn't exist.")
            return
        else:
            logging.info("Loading data")
            ##      Data loading
            if self.source.endswith('.csv'):
                sep = ','
            elif self.source.endswith('.tsv'):
                sep = '\t'
            else:
                raise ValueError('Unsupported file format')

            self._dataframe = pd.read_csv(self.source, sep=sep)

    @ray_manager()
    def __standardize(self):
        """
        This method runs the standardization of the molecules represented in a dataframe. It requires to have a molecular representation (either
        the smiles string or the RDKit mol object). If the data has already been standardized, does nothing (thus its safe to call it at any point).
        """
        if self.is_standardized:
            return

        self._dataframe = self._dataframe.reset_index(drop=True)
        self._dataframe = self._dataframe.rename(columns={self._dataframe.columns[self._dataframe.columns.str.lower().tolist().index(self.NAME_SMILES_RAW.lower())]: self.NAME_SMILES_RAW})

        n_processes = min(N_CPUS, os.cpu_count())  # TODO: improve parallelization (generalize)
        self._dataframe[self.NAME_MOLOBJ] = parallel_apply(
            dataframe=self._dataframe[self.NAME_SMILES_RAW],
            function=molFromSmiles,
            n_processes=n_processes,
            desc="Building mol objects",
        )
        self._dataframe[self.NAME_MOLOBJ] = parallel_apply(
            dataframe=self._dataframe[self.NAME_MOLOBJ],
            function=lambda mol: standardize(mol) if not mol is None else None,
            n_processes=n_processes,
            desc="Standardizing molecules",
        )
        self._dataframe[self.NAME_INCHIKEY] = parallel_apply(
            dataframe=self._dataframe[self.NAME_MOLOBJ],
            function=lambda mol: Chem.MolToInchiKey(mol) if not mol is None else None,
            n_processes=n_processes,
            desc="Retrieving InChiKey",
        )
        self._dataframe[self.NAME_SMILES_STD] = parallel_apply(
            dataframe=self._dataframe[self.NAME_MOLOBJ],
            function=lambda mol: Chem.MolToSmiles(mol) if not mol is None else None,
            n_processes=n_processes,
            desc="Retrieving SMILES",
        )

        self.is_standardized = True

    def _deduplicate(self, endpoint2task:dict={}, subset:list=[]):
        """
        This method identifies whether there are duplicated molecules in the dataframe, combining the molecular property of the same molecules
        into a single value.

        Standardization must have run successfully first, so that needed columns (INCHIKEY, SMILES, MOLOBJ) are available.

        !It cannot be reverted!
        """
        if self.is_deduplicated:
            return

        endpoint2task = {endpoint:task.upper() for endpoint,task in endpoint2task.items()}  # Ensure that the provided tasks are uppercase

        if len(subset) == 0:
            subset.append(self.NAME_INCHIKEY)

        ##      Ensure standardization
        self.__standardize()

        original_size = len(self._dataframe)
        if self._dataframe.duplicated(subset).any():
            logging.info("Deduplicating multiple value annotations for a single molecule-endpoint")
            deduplicate_funcs = {  # TODO: change name
                self.NAME_ID : (lambda x: tuple(x)),  # TODO: there could be multiple ids here
                self.NAME_MOLOBJ : (lambda x: list(set(x))[0]),
                self.NAME_SMILES_RAW : (lambda x: tuple(x)),  # TODO: there could be multiple smiles here
                self.NAME_SMILES_STD : (lambda x: list(set(x))[0]),
                self.NAME_REF : (lambda x: ", ".join(set(x).difference({np.nan}))),
            }
            # all columns that are in the original dataframe and not have a specific funcion assigned for deduplication are saved as a tuple
            def unique_tuple(x):
                try:
                    return tuple(set(x))
                except:
                    return tuple(x)
            # deduplicate_funcs.update({col:(lambda x: tuple(x)) for col in self._dataframe.columns if col not in deduplicate_funcs.keys()})
            deduplicate_funcs.update({col:(lambda x: unique_tuple(x)) for col in self._dataframe.columns if col not in deduplicate_funcs.keys()})

            if self.NAME_ENDPOINT in self._dataframe.columns:
                groups = []
                for endpoint, group in self._dataframe.groupby(self.NAME_ENDPOINT):
                    funcs2apply = {col:func for col, func in deduplicate_funcs.items() if col in set(group.columns).difference(subset)}
                    # always defaults to regression
                    if self.NAME_VALUE in self._dataframe.columns:
                        if endpoint2task.get(endpoint, TASK_REGRESSION) == TASK_CLASSIFICATION:
                            funcs2apply.update({self.NAME_VALUE: (lambda x: list(set(x))[0] if len(set(x))==1 else None)})
                        else:
                            funcs2apply.update({self.NAME_VALUE: (lambda x: x.mean() if len(x) < 5 else x[x.between(x.mean() - 2*x.std(), x.mean() + 2*x.std())].mean())})
                    group = group.groupby(list(set(subset).intersection(group.columns))).agg(funcs2apply).reset_index().dropna(subset=list(funcs2apply.keys()))
                    group[self.NAME_ENDPOINT] = endpoint
                    groups.append(group)
                self._dataframe = pd.concat(groups).reset_index(drop=True)
            else:   # if no values/endpoints are specified
                funcs2apply = {col:func for col, func in deduplicate_funcs.items() if col in set(self._dataframe.columns).difference(subset)}
                self._dataframe = self._dataframe.groupby(list(set(subset).intersection(self._dataframe.columns))).agg(funcs2apply).reset_index().dropna(subset=list(funcs2apply.keys()))
        else:  # no duplicates in subset
            for col in [self.NAME_ID, self.NAME_SMILES_RAW]:
                if col in self._dataframe.columns:
                    self._dataframe[col] = self._dataframe[col].apply(lambda x: tuple([x]))
        # if some columns are converted to single-item tuples, explode them
        single_item_columns = [col for col in self._dataframe.columns if all(self._dataframe[col].apply(lambda x: len(x) if isinstance(x, tuple) else np.nan)==1)]
        if single_item_columns:
            self._dataframe = self._dataframe.explode(column=single_item_columns)

        if len(self._dataframe) < original_size:
            logging.info(f"Reduced to {len(self._dataframe)} single molecule-endpoints (from {original_size})")
            self.__buildMolecules()
        
        self.is_deduplicated = True
    
    def __buildMolecules(self):
        """
        Builds the molecules from the dataframe, if the second is available.
        """

        if self._dataframe is None:
            logging.warning("Trying to build molecules from empty DF")
            return
        
        self._molecules = {}

        logging.info("Building Molecules")
        for _, row in tqdm(self._dataframe.iterrows(), total=len(self._dataframe), desc="Adding molecules"):
            molObj = row[self.NAME_MOLOBJ] if self.NAME_MOLOBJ in row else Chem.MolFromSmiles(row[self.NAME_SMILES_STD])
            self.__addMolecule(
                id=row[self.NAME_ID],
                smiles=row[self.NAME_SMILES_STD],
                inchikey=row[self.NAME_INCHIKEY],
                molObj=molObj,
                standard=self.is_standardized,
                # Add other fields as necessary,
            )

        logging.info(f"Molecules ready ({len(self._molecules)} molecules)")

    def __addMolecule(self, **props):
        if self._molecules is None:
            self._molecules = {}
        molecule_info = MoleculeInfo(**props)
        if isinstance(props[self.NAME_ID], tuple):
            self._molecules[props[self.NAME_ID]] = molecule_info
        else:
            self._molecules[tuple([props[self.NAME_ID]])] = molecule_info
    
    def Molecules(self):

        if not self._molecules is None:
            return self._molecules
        
        ##  Otherwise, build (if DF is available)
        if not self._dataframe is None:
            self.__buildMolecules()
            return self._molecules
    
    def __calcFeatureDF(self, feature):
        """
        This is one of the AVAILABLE FEATURE EXTRACTORS for molecules.

        Once this level of depth in the code is reached, features will be COMPUTED.

        Runs over the molecules in list, extracting features and builds a DF with one molecule per row and each feature represented 
        in a column of the dataframe. The primary key of this dataframe is the molecule identifier.
        """

        logging.info(f"Building {feature}")

        molecules = self.Molecules()

        ##      feature as a list of lists
        if feature in self.featuresInDF_list:
            feature_list = [[id]+list(molecule.getFeature(feature).values()) for id, molecule in molecules.items()]
        else:
            feature_list = parallel_apply(pd.Series(molecules), function=(lambda mol: [mol.id]+list(mol.getFeature(feature).values())), desc=f"Calculating {feature}").tolist()

        feature_df = pd.DataFrame(feature_list, columns=[MoleculeData.NAME_ID]+[f"{feature}_{feat}" for feat in MoleculeInfo.AVAILABLE_FEATURES[feature]])

        return feature_df
    
    def __addFeatureToDF(self, feature):
        """
        This method checks whether the feature has been added to the dataframe (via the featuresInDF list), and either computes it
        or just returns the list of columns with their names, for higher levels to handle the dataframe properly.

        At this level, the feature may or may not be computed, depending on whether it has been calculated before hand.
        """
        if isinstance(feature, str):
            feature = feature.upper()
            if not feature in MoleculeInfo.AVAILABLE_FEATURES:
                logging.error(f"Feature {feature} not available, choose one from {list(MoleculeInfo.AVAILABLE_FEATURES.keys())}")
                return
            if not feature in self.featuresInDF_list:
                feature_df = self.__calcFeatureDF(feature).round(VALUES_ROUNDING)
                self._dataframe = self._dataframe.merge(feature_df, on=self.NAME_ID)
            feature_columns = [f"{feature}_{col}" for col in MoleculeInfo.AVAILABLE_FEATURES[feature]]
        elif isinstance(feature, pd.DataFrame):
            feature_df = feature.copy()
            feature = "CUSTOM"
            feature_columns = [f"{feature}_{col}" if col != self.NAME_ID else self.NAME_ID for col in feature_df.columns]
            feature_df.columns = feature_columns
            if self.NAME_ID in feature_columns:
                feature_columns.remove(self.NAME_ID)
            else:
                feature_df[self.NAME_ID] = self._dataframe[self.NAME_ID].explode().values
            feature_df = feature_df.set_index(self.NAME_ID).loc[self._dataframe[self.NAME_ID].apply(lambda l: l[0]).values].set_index(self._dataframe[self.NAME_ID]).reset_index()
            self._dataframe = self._dataframe.merge(feature_df, on=self.NAME_ID, suffixes=("_to_drop", ""), how="left")
            self._dataframe = self._dataframe.drop(columns=[col for col in self._dataframe if col.endswith("_to_drop")])
        else:
            logging.error(f"Feature {feature} format not recognized")

        self.featuresInDF_list.append(feature)

        return feature_columns
    
    def DataFrame(self, features:list = [], columns:list = []):

        '''
        Overloaded method from DataForModeling. Returns the dataframe containing the desired molecule data, filtering out by the required columns
        and features. Features must be defined according to specific names provided by this class (see above).
        '''
        # configsObj = Configs()
        
        ##      If no features are given, return DF as is
        if (not self._dataframe is None) and (len(features) == 0) and (len(columns) == 0):
            return self._dataframe
        ##      Delete any specific requirement that is not reflected in the dataframe
        columns = list(set(columns).intersection(self._dataframe.columns.tolist()))
        # features_columns = []

        ##      Extract the required features, using the logic flag to avoid reconstructing them (TODO)
        return_df = self._dataframe   #   Default: Return full dataframe

        ##      Add Features
        for feature in features:
            feature_cols = self.__addFeatureToDF(feature)
            columns.extend(feature_cols)
        
        ##      Select columns based on the given argument
        if (len(columns) > 0):
            return_df = self._dataframe[columns]

        return return_df
    
    ###
    def splitBy(self, column):
        return {col:self.copy(_dataframe=group) for col, group in self.DataFrame().groupby(column)}

    def copy(self, **props):
        """
        Create a copy of the current `MoleculeData` object, optionally modifying specific attributes.

        Parameters:
        ----------
        props : keyword arguments, optional
            Keyword arguments representing property names and values to override in the new object. 
            If no keyword arguments are provided, the original object's properties are retained.

        Returns:
        -------
        MoleculeData
            A new `MoleculeData` object with copied properties and molecules.

        Behavior:
        --------
        - All attributes of the original object (`self`) are deep-copied into the new object unless overridden by 
        the provided keyword arguments (`props`).
        - Molecules from the original object are mapped to the new object's dataframe based on the `NAME_ID` attribute.
          Any molecules present in the dataframe but missing from the original `molecules` dictionary trigger the creation
          of new molecule entries. Molecules in the original `molecules` dictionary but missing from the dataframe are removed.

        Example:
        -------
        new_mol = mol.copy(is_standardized=True)
        """

        newObj = MoleculeData(source=None)

        # Loop through all properties in the current object's dictionary.
        # If the property is in `props`, override the value with the one provided.
        for prop, value in self.__dict__.items():
            if prop in props:
                value = deepcopy(props[prop])
            if hasattr(newObj, prop):
                setattr(newObj, prop, deepcopy(value))

        if self._molecules is None:
            newObj._molecules = None
        else:
            newObj._molecules = {id:self._molecules.get(id) for id in newObj._dataframe[self.NAME_ID]}
            missing_mols = newObj._dataframe[~newObj._dataframe[self.NAME_ID].isin(newObj._molecules.keys())]
            # If there are molecules that exist in the new dataframe but are missing in the copied `molecules`, add them to the new `molecules``.
            if len(missing_mols):
                logging.info("New molecules in dataframe, Building molecules")
                for _, row in tqdm(missing_mols.iterrows(), total=len(missing_mols)):
                    molObj = row[self.NAME_MOLOBJ] if self.NAME_MOLOBJ in row else Chem.MolFromSmiles(row[self.NAME_SMILES_STD])
                    newObj.__addMolecule(
                        id=row[self.NAME_ID],
                        smiles=row[self.NAME_SMILES_STD], 
                        inchikey=row[self.NAME_INCHIKEY],
                        molObj=molObj, 
                        standard=self.is_standardized,
                        # Add other fields as necessary,
                    )
            # Update item ids
            try:
                # if dataframe deduplicated ids are tuples
                newObj._ids = [id for ids in newObj._molecules.keys() for id in ids]
            except:
                newObj._ids = list(newObj._molecules.keys())

        return newObj

    ###
    def __len__(self):
        return len(self.DataFrame())