import pytest
import pandas as pd
import numpy as np
from scipy.stats import skewtest
from rdkit import Chem
from assay_inspector.AI_Calculation import Calculation
from assay_inspector.AI_MoleculeData import MoleculeData
from assay_inspector.AI_MoleculeInfo import MoleculeInfo

# --- Mock Classes and Fixtures for Testing ---
class MockMainCalculation:
    def __init__(self, task="REGRESSION", feature_type="rdkit", outliers_method="zscore", distance_metric="euclidean"):
        self.endpoint_name = "test_endpoint"
        self.task = task
        self.feature_type = feature_type
        self.outliers_method = outliers_method
        self.distance_metric = distance_metric
        self.lower_bound = -10
        self.upper_bound = 10
        self.sources_list = ['source1', 'source2']
        self.hex_colors = ['#1f77b4', '#ff7f0e']

        if self.feature_type == 'rdkit':
            self.features = MoleculeInfo.FEAT_RDKIT_DESC
        elif self.feature_type == 'ecfp4':
            self.features = MoleculeInfo.FEAT_ECFP4
        else: # for 'custom' and other cases
            self.features = "CUSTOM"

@pytest.fixture
def mock_calculation():
    return Calculation(mainSelf=MockMainCalculation())

# --- Skewness Calculation Tests ---
def test_calculate_skewness_right_skewed(mock_calculation):
    """Test that calculateSkewness correctly identifies a right-skewed distribution."""
    np.random.seed(42)
    right_skewed_data = pd.DataFrame({
        'value': np.random.lognormal(mean=0, sigma=1, size=10000)
    })

    skewness, statistic, pvalue = mock_calculation.calculateSkewness(right_skewed_data)
    assert skewness > 0
    assert pvalue < 0.05
    assert isinstance(skewness, float)
    assert isinstance(statistic, float)
    assert isinstance(pvalue, float)

def test_calculate_skewness_normal(mock_calculation):
    """Test that calculateSkewness correctly identifies a normal distribution."""
    np.random.seed(42)
    normal_data = pd.DataFrame({
        'value': np.random.normal(loc=0, scale=1, size=10000)
    })

    skewness, statistic, pvalue = mock_calculation.calculateSkewness(normal_data)
    assert abs(skewness) < 0.1
    assert pvalue > 0.05

# --- Kurtosis Calculation Test ---
def test_calculate_kurtosis_normal():
    """Test that calculateKurtosis correctly identifies a normal distribution."""
    mock_calc = Calculation(mainSelf=MockMainCalculation())
    np.random.seed(42)
    normal_data = pd.DataFrame({
        'value': np.random.normal(loc=0, scale=1, size=10000)
    })

    kurtosis_val, statistic, pvalue = mock_calc.calculateKurtosis(normal_data)
    assert pytest.approx(kurtosis_val, abs=0.1) == 0.0
    assert pvalue > 0.05

# --- Outlier Identification (IQR method) Test ---
def test_identify_outliers_iqr():
    """Test that identifyOutliers correctly finds outliers using the IQR method."""
    mock_calc_iqr = Calculation(mainSelf=MockMainCalculation(outliers_method="iqr"))

    data = {'value': np.concatenate((np.random.normal(loc=10, scale=2, size=100), [200]))}
    data_df = pd.DataFrame(data)
    data_df['inchikey'] = [f'mol_{i}' for i in range(len(data_df))]
    data_df['endpoint'] = "test"
    data_df['ref'] = "source1"

    outliers_info, outliers_set = mock_calc_iqr.identifyOutliers(data_df)
    assert outliers_info[4] == 1
    assert outliers_set == {'mol_100'}

# --- Distance Matrix Calculation Test ---
def test_calculate_distance_matrix_ecfp4():
    """Test that calculateDistanceMatrix correctly generates a Tanimoto distance matrix for ECFP4."""
    # Use molecules that share a common scaffold to ensure a non-zero similarity
    mol_data = MoleculeData(source=pd.DataFrame({'smiles': ['c1ccccc1', 'Cc1ccccc1']}))
    features_df = mol_data.DataFrame(features=[MoleculeInfo.FEAT_ECFP4])
    
    mock_calc_ecfp4 = Calculation(mainSelf=MockMainCalculation(feature_type="ecfp4"))
    distance_matrix = mock_calc_ecfp4.calculateDistanceMatrix(features_df)

    assert distance_matrix.shape == (2, 2)
    assert distance_matrix[0, 1] > 0