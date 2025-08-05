import pytest
import pandas as pd
import numpy as np
from assay_inspector.AI_FeaturesPreprocessing import FeaturesPreprocessing
from assay_inspector.AI_MoleculeData import MoleculeData

# --- Mock Classes for Testing ---
# This mock class provides the necessary attributes for the FeaturesPreprocessing class to work.
class MockMainFeatures:
    def __init__(self):
        self.endpoint_name = "test_endpoint"
        self.task = "REGRESSION"
        self.feature_type = "rdkit"
        self.features = 'RDKIT'

@pytest.fixture
def mock_preprocessor():
    return FeaturesPreprocessing()

# --- Test Case 1: Imputation of Missing Values ---
def test_impute_missing_values(mock_preprocessor):
    """Test that FeaturesPreprocessing correctly imputes missing values with the mean strategy."""
    # Start with clean data
    data = {'smiles': ['CCO', 'CCC'], 'value': [1.0, 2.0]}
    df = pd.DataFrame(data)

    mol_data = MoleculeData(source=df)
    
    # Let MoleculeData compute the features, then introduce the NaN value
    full_df = mol_data.DataFrame(features=['RDKIT'], columns=['id', 'inchikey'])
    mol_data = mol_data.copy(_dataframe=full_df)
    mol_data.DataFrame().loc[1, 'RDKIT_MolWt'] = np.nan
    
    preprocessed_data = mock_preprocessor.fit_transform(mol_data, features_ids=['RDKIT'], endpoint2task={'test': 'REGRESSION'})
    preprocessed_df = preprocessed_data.DataFrame()
    
    # The imputed value is the mean (44.08). The min value is 44.08. So the scaled value should be 0.0
    assert not preprocessed_df['RDKIT_MolWt'].isnull().any()
    assert pytest.approx(preprocessed_df['RDKIT_MolWt'].iloc[1]) == 0.0

# --- Test Case 2: Feature Scaling (MinMaxScaler) ---
def test_scale_features(mock_preprocessor):
    """Test that FeaturesPreprocessing correctly scales features using MinMaxScaler."""
    # Start with clean data
    data = {'smiles': ['CCO', 'CCC'], 'value': [1.0, 2.0]}
    df = pd.DataFrame(data)

    mol_data = MoleculeData(source=df)
    
    # Let MoleculeData compute the features
    full_df = mol_data.DataFrame(features=['RDKIT'], columns=['id', 'inchikey'])
    mol_data = mol_data.copy(_dataframe=full_df)

    preprocessed_data = mock_preprocessor.fit_transform(mol_data, features_ids=['RDKIT'], endpoint2task={'test': 'REGRESSION'})
    preprocessed_df = preprocessed_data.DataFrame()
    
    # With MinMaxScaler, the values should be scaled between 0 and 1
    assert pytest.approx(preprocessed_df['RDKIT_MolWt'].min()) == 0.0
    assert pytest.approx(preprocessed_df['RDKIT_MolWt'].max()) == 1.0