import pytest
from assay_inspector.AI_MoleculeInfo import MoleculeInfo
from rdkit import Chem

# --- Test Case 1: RDKit Descriptors ---
def test_get_rdkit_descriptors():
    """Test that getFeature correctly computes RDKit descriptors for a molecule."""
    smiles = "c1ccccc1" # Benzene
    mol_info = MoleculeInfo(smiles=smiles)
    descriptors = mol_info.getFeature(MoleculeInfo.FEAT_RDKIT_DESC)

    # Check that the number of descriptors matches the available ones
    assert len(descriptors) == len(MoleculeInfo.AVAILABLE_FEATURES[MoleculeInfo.FEAT_RDKIT_DESC])
    
    # Check that a known descriptor is present and has a valid value
    assert "MolWt" in descriptors
    assert pytest.approx(descriptors["MolWt"], 0.01) == 78.114

# --- Test Case 2: ECFP4 Fingerprints ---
def test_get_ecfp4_fingerprints():
    """Test that getFeature correctly computes ECFP4 fingerprints."""
    smiles = "CCO" # Ethanol
    mol_info = MoleculeInfo(smiles=smiles)
    fingerprints = mol_info.getFeature(MoleculeInfo.FEAT_ECFP4)

    # Check that the number of bits matches the hardcoded value
    assert len(fingerprints) == MoleculeInfo.ECFP4_nBits
    
    # Check that the bits are binary (0 or 1)
    assert all(bit in [0, 1] for bit in fingerprints.values())