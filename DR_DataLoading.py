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

###
def __getMoleculeProfileDataframe(self, data, reporting):

    """
    Generates the DataFrame containing the feature profiles of the molecules. 
    """

    # Deduplicate the dataset 
    if reporting == 'individual':
        data._deduplicate(subset=[CONFIGS.NAMES.INCHIKEY], endpoint2task={self.endpoint_name:self.task})
    elif reporting == 'comparative':
        data._deduplicate(subset=[CONFIGS.NAMES.INCHIKEY,CONFIGS.NAMES.REF], endpoint2task={self.endpoint_name:self.task})

    # Select the data corresponding to the given endpoint
    data = data.splitBy(CONFIGS.NAMES.ENDPOINT_ID)[self.endpoint_name]

    # Perform feature preprocessing
    preprocessing = FeaturesPreprocessing()
    data = preprocessing.fit_transform(data, features_ids=[self.features], endpoint2task={self.endpoint_name:self.task})

    # Get the molecule DataFrame
    data_df = data.getDataframe(features=[self.features],
                                columns=[CONFIGS.NAMES.INCHIKEY, CONFIGS.NAMES.VALUE, CONFIGS.NAMES.REF, CONFIGS.NAMES.ENDPOINT_ID])

    return data_df

###
def __getInChIKeySet(self, mol_data):

    """
    Generates the set of InChIKeys set for the molecules in the reference set. 
    """

    # Get the DataFrame of reference molecules
    ref_df = mol_data.getDataframe(columns=[CONFIGS.NAMES.INCHIKEY, 'smiles']) # TODO: Change to a constant?
    
    # Extract the InChIKey set
    inchikeys_set = set(ref_df[CONFIGS.NAMES.INCHIKEY].tolist())

    return inchikeys_set