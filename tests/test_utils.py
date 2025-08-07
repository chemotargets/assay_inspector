import pytest
from assay_inspector import AI_Utils
from rdkit import Chem

def test_mol_from_smiles_valid():
    """Test that molFromSmiles correctly creates a molecule object for a valid SMILES string."""
    smiles = "c1ccccc1"
    mol = AI_Utils.molFromSmiles(smiles)
    assert isinstance(mol, Chem.Mol)
    assert Chem.MolToSmiles(mol) == smiles

def test_mol_from_smiles_invalid():
    """Test that molFromSmiles returns None for an invalid SMILES string."""
    smiles = "InvalidSmiles"
    mol = AI_Utils.molFromSmiles(smiles)
    assert mol is None

def test_standardize_removes_salt():
    """Test that standardize correctly removes a salt counter-ion."""
    # A molecule with a salt counter-ion
    smiles_with_salt = "C(Cl)C(Cl)C(=O)[O-].[Na+]"
    mol_with_salt = Chem.MolFromSmiles(smiles_with_salt)
    standardized_mol = AI_Utils.standardize(mol_with_salt)
    standardized_smiles = Chem.MolToSmiles(standardized_mol)

    # The expected standardized smiles should only contain the molecule, not the salt
    assert standardized_smiles == "O=C([O-])C(Cl)CCl"