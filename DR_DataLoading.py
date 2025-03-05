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

from DR_FeaturesPreprocessing import FeaturesPreprocessing

import json

### Configs
with open("configs.json", "r") as file:
    config = json.load(file)

###
def __getMoleculeProfileDataframe(self, data, reporting):

    """
    Generates the DataFrame containing the feature profiles of the molecules. 
    """

    # Deduplicate the dataset 
    if reporting == 'individual':
        data._deduplicate(subset=[config["NAMES"]["INCHIKEY"]], endpoint2task={self.endpoint_name:self.task})
    elif reporting == 'comparative':
        data._deduplicate(subset=[config["NAMES"]["INCHIKEY"],config["NAMES"]["REF"]], endpoint2task={self.endpoint_name:self.task})

    # Select the data corresponding to the given endpoint
    data = data.splitBy(config["NAMES"]["ENDPOINT_ID"])[self.endpoint_name]

    # Perform feature preprocessing
    preprocessing = FeaturesPreprocessing()
    data = preprocessing.fit_transform(data, features_ids=[self.features], endpoint2task={self.endpoint_name:self.task})

    # Get the molecule DataFrame
    data_df = data.getDataframe(features=[self.features],
                                columns=[config["NAMES"]["INCHIKEY"], config["NAMES"]["VALUE"], config["NAMES"]["REF"], config["NAMES"]["ENDPOINT_ID"]])

    return data_df

###
def __getInChIKeySet(self, mol_data):

    """
    Generates the set of InChIKeys set for the molecules in the reference set. 
    """

    # Get the DataFrame of reference molecules
    ref_df = mol_data.getDataframe(columns=[config["NAMES"]["INCHIKEY"], config["NAMES"]["SMILES"]]) 
    
    # Extract the InChIKey set
    inchikeys_set = set(ref_df[config["NAMES"]["INCHIKEY"]].tolist())

    return inchikeys_set