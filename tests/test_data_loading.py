import pytest
import pandas as pd
from assay_inspector.AI_MoleculeData import MoleculeData

def test_deduplicate_by_inchikey_and_ref():
    """Test that deduplication correctly removes duplicates based on inchikey and ref."""
    # The correct InChIKey for the SMILES 'CCO'
    inchikey_cco = 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'

    # Create a DataFrame with some duplicate and unique entries
    data = [
        {'id': 1, 'smiles': 'CCO', 'inchikey': 'XYZ', 'ref': 'source1', 'value': 1.0, 'endpoint': 'test'},
        {'id': 2, 'smiles': 'CCO', 'inchikey': 'XYZ', 'ref': 'source1', 'value': 2.0, 'endpoint': 'test'}, # Duplicate of the first
        {'id': 3, 'smiles': 'CCO', 'inchikey': 'XYZ', 'ref': 'source2', 'value': 3.0, 'endpoint': 'test'}, # Unique, as ref is different
        {'id': 4, 'smiles': 'CC', 'inchikey': 'ABC', 'ref': 'source1', 'value': 4.0, 'endpoint': 'test'}, # Unique
    ]
    df = pd.DataFrame(data)

    # Instantiate the MoleculeData class and run deduplication
    mol_data = MoleculeData(source=df)
    mol_data._deduplicate(subset=[MoleculeData.NAME_INCHIKEY, MoleculeData.NAME_REF], endpoint2task={'test': 'REGRESSION'})
    deduplicated_df = mol_data.DataFrame()

    # The expected number of rows after deduplication is 3
    assert len(deduplicated_df) == 3
    
    # Assert that the merged row for 'source1' has an average value
    # For regression, duplicates are averaged (1.0 + 2.0) / 2 = 1.5
    assert round(deduplicated_df.loc[(deduplicated_df['inchikey'] == inchikey_cco) & (deduplicated_df['ref'] == 'source1'), 'value'].iloc[0], 2) == 1.5

def test_deduplicate_without_endpoint():
    """Test deduplication when there are no endpoint values defined."""
    # The correct InChIKey for the SMILES 'CCO'
    inchikey_cco = 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'

    data = [
        {'id': 1, 'smiles': 'CCO', 'inchikey': 'XYZ', 'ref': 'source1'},
        {'id': 2, 'smiles': 'CCO', 'inchikey': 'XYZ', 'ref': 'source1'}, # Duplicate
        {'id': 3, 'smiles': 'CCO', 'inchikey': 'XYZ', 'ref': 'source2'}, # Unique
    ]
    df = pd.DataFrame(data)
    mol_data = MoleculeData(source=df)
    mol_data._deduplicate(subset=[MoleculeData.NAME_INCHIKEY, MoleculeData.NAME_REF])
    deduplicated_df = mol_data.DataFrame()
    
    assert len(deduplicated_df) == 2
    assert inchikey_cco in deduplicated_df['inchikey'].values

def test_load_dataframe_correctly():
    """Test that the MoleculeData class correctly loads a pandas DataFrame."""
    data = {'smiles': ['CCO', 'CCC'], 'value': [1.0, 2.0]}
    df = pd.DataFrame(data)

    mol_data = MoleculeData(source=df)
    
    assert mol_data.DataFrame() is not None
    assert len(mol_data.DataFrame()) == 2
    assert 'molObj' in mol_data.DataFrame().columns
    assert mol_data.is_standardized == True

