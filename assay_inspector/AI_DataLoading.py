#!/usr/bin/env python

__author__ = "Raquel Parrondo-Pizarro"
__date__ = "20250220"
__copyright__ = "Copyright 2024, Chemotargets"
__license__ = ""
__credits__ = ["Data Science & Translational Research Group"]
__maintainer__ = "Raquel Parrondo-Pizarro"
__version__ = "20250220"
__deprecated__ = False

### Imports

from AI_FeaturesPreprocessing import FeaturesPreprocessing

### Configs
import os
import json

config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'configs.json')

with open(config_path, "r") as file:
    config = json.load(file)

class DataLoading():
    """
    Class to load and generate the molecule dataframe for both endpoint dataset and reference
    set (if provided). 
    """

    def __init__(self, mainSelf):
        self.__mainSelf = mainSelf

    ###
    def getMoleculeProfileDataframe(self, data, reporting):

        """
        Generates the DataFrame containing the feature profiles of the molecules. 
        """        
        # Deduplicate the dataset 
        if reporting == 'individual':
            data._deduplicate(subset=[config["NAMES"]["INCHIKEY"]], endpoint2task={self.__mainSelf.endpoint_name:self.__mainSelf.task})
        elif reporting == 'comparative':
            data._deduplicate(subset=[config["NAMES"]["INCHIKEY"],config["NAMES"]["REF"]], endpoint2task={self.__mainSelf.endpoint_name:self.__mainSelf.task})

        # Select the data corresponding to the given endpoint
        data = data.splitBy(config["NAMES"]["ENDPOINT_ID"])[self.__mainSelf.endpoint_name]

        if self.__mainSelf.feature_type == 'custom':
            if not data.NAME_ID in self.__mainSelf.features.columns:
                self.__mainSelf.features[data.NAME_ID] = data.DataFrame()[data.NAME_ID]

        # Perform feature preprocessing
        preprocessing = FeaturesPreprocessing()
        data = preprocessing.fit_transform(data, features_ids=[self.__mainSelf.features], endpoint2task={self.__mainSelf.endpoint_name:self.__mainSelf.task})

        # Get the molecule DataFrame
        # data_df = data.DataFrame(features=data.featuresInDF_list,
        #                          columns=[config["NAMES"]["INCHIKEY"], config["NAMES"]["VALUE"], config["NAMES"]["REF"], config["NAMES"]["ENDPOINT_ID"]])
        data_df = data.DataFrame()[[config["NAMES"]["INCHIKEY"], config["NAMES"]["VALUE"], config["NAMES"]["REF"], config["NAMES"]["ENDPOINT_ID"]]+[col for col in data.DataFrame().columns for feat in data.featuresInDF_list if col.startswith(feat)]]

        return data_df

    ###
    def getInChIKeySet(self, mol_data):

        """
        Generates the set of InChIKeys set for the molecules in the reference set. 
        """

        # Get the DataFrame of reference molecules
        ref_df = mol_data.DataFrame(columns=[config["NAMES"]["INCHIKEY"], config["NAMES"]["SMILES"]]) 
        
        # Extract the InChIKey set
        inchikeys_set = set(ref_df[config["NAMES"]["INCHIKEY"]].tolist())

        return inchikeys_set