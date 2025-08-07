#! /usr/bin/env python

__author__ = "Luca Menestrina"
__date__ = "20250130"
__copyright__ = "Copyright 2025, Chemotargets"
__license__ = ""
__credits__ = ["Data Science & Translational Research Group", "Ismael Gomez", "Ricard Garcia"]
__maintainer__ = "Luca Menestrina"
__version__ = "20250224"
__deprecated__ = False

from rdkit import Chem

from .AI_Utils import logging, standardize
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
# from mordred import Calculator as mordred_calculator, descriptors as mordred_descriptors

class MoleculeInfo():

    """
    This class handles all possible molecular representations that may be needed through the PA engine execution process.
    """

    ##      Feature names - These do not correspond to actual DF columns, but with features that are represented by many columns within a DF.
    FEAT_ECFP4 = "ECFP4"
    ECFP4_nBits = 1024  # TODO: should it be hardcoded?
    ECFP4_radius = 2  # TODO: should it be hardcoded?
    FEAT_RDKIT_DESC = "RDKIT"
    # FEAT_MORDRED_DESC = "MORDRED"
    FEAT_CUSTOM = "CUSTOM"

    AVAILABLE_FEATURES = {
        FEAT_ECFP4 : list(range(1,ECFP4_nBits+1)),
        FEAT_RDKIT_DESC : [x[0] for x in Descriptors._descList if x[0] not in ["Ipc"]],
        # FEAT_MORDRED_DESC : [str(desc) for desc in mordred_calculator(mordred_descriptors, ignore_3D=False).descriptors],
    }

    ###
    def __init__(self, id = None, smiles = None, inchikey = None, molObj = None, standard = False):
    
        """
        Class constructor, initializes the whole set of molecular repesentations (empty at first).

        Arguments:
        -   id:         Custom molecule identifier
        -   smiles:     String with Smiles representation of the molecule.
        -   inchikey:   String with Inchi key of the molecule.
        -   molObj:     RDKit mol object.
        -   standard:   Boolean representing if the molecule has already been standardized.
        """

        self.id = id
        self.smiles = smiles
        self.inchikey = inchikey
        self.molObj = molObj

        ##      Flag for standardization
        self.standard = standard

        self.features = {}

        # ##      RDKit descriptors (format: dict of descriptor_id: descriptor_val)
        # self.descriptorsRDKit = None

        ##      Standardization (must be done when possible)
        if not smiles is None:
            self.molObj = Chem.MolFromSmiles(self.smiles)

        if (not self.molObj is None) and (not standard):
            self.__Standardize()                                # In theory should work by standardizing the object internally
            if not self.molObj is None:
                self.smiles = Chem.MolToSmiles(self.molObj)         # Updated smiles
                self.inchikey = Chem.MolToInchiKey(self.molObj)     # Updated inchikey

    ###
    def __Standardize(self):
        """
        Standardizes the mol object if available. Sets up the standard flag to True, so this is known.
        """
        
        if (self.standard) or (self.molObj is None):
            return
        
        self.molObj = standardize(self.molObj)      # TODO: Use setter?
        self.standard = True                        # TODO: Handle based on return


    ###
    @property
    def getID(self):

        return self.id

    ###
    @property
    def getInchiKey(self):

        """
        TODO: Compute inchikey if not available
        """
        return self.inchikey
    
    ###
    @property
    def mol(self):
        
        """
        Returns mol object (computes it if missing)
        """

        if not self.molObj is None:
            self.__Standardize()                            
            return self.molObj
        
        ##  Mol obj computing: (1) from smiles
        if not self.smiles is None:
            self.molObj = Chem.MolFromSmiles(self.smiles)
            self.__Standardize()
            return self.molObj
        
    ###
    def _calcFeature(self, feature):
        """
        Utility function for computing molecular features.

        This function calculates molecular features for a given molecule.

        Parameters:
            - feature: the feature group name to compute (RDKIT, ECFP4, ...)

        Returns:
            - a dictionary mapping feature names (individual) to their calculated values
        """
        
        if not feature in self.AVAILABLE_FEATURES:
            logging.error(f"Feature {feature} not available, choose one from {list(self.AVAILABLE_FEATURES.keys())}")
            return
        if self.molObj is None:
            feature_vals = [None]*len(self.AVAILABLE_FEATURES[feature])
        else:
            if feature == self.FEAT_ECFP4:
                feature_vals = AllChem.GetMorganFingerprintAsBitVect(self.mol, radius=self.ECFP4_radius, nBits=self.ECFP4_nBits).ToList()
            elif feature == self.FEAT_RDKIT_DESC:
                # TODO: filter self.AVAILABLE_FEATURES[feature] checking if others have already computed them
                calc = MoleculeDescriptors.MolecularDescriptorCalculator(self.AVAILABLE_FEATURES[feature])
                feature_vals = calc.CalcDescriptors(self.mol)
            # elif feature == self.FEAT_MORDRED_DESC:
            #     calc = mordred_calculator(mordred_descriptors, ignore_3D=False)(self.mol)
            #     feature_vals = calc.values()
            else:
                logging.error(f"Feature {feature} not available, choose one from {list(self.features.keys())}")

        feature_dict = dict(zip(self.AVAILABLE_FEATURES[feature], feature_vals))

        self.features[feature] = feature_dict

        return feature_dict

    ###
    def getFeature(self, feature):
        if not feature in self.AVAILABLE_FEATURES:
            logging.error(f"Feature {feature} not available, choose one from {list(self.AVAILABLE_FEATURES.keys())}")
            return

        if not self.features.get(feature, None) is None:
            return self.features[feature]
        
        return self._calcFeature(feature)