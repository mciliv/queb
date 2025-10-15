import pytest
from pathlib import Path
from rdkit import Chem
from sdf import sdf

@pytest.fixture
def cleanup_sdf():
    """Fixture to clean up any generated SDF files after the test."""
    yield
    for sdf_file in Path("tests/sdf_files").glob("*.sdf"):
        if sdf_file.name.startswith("test_"):
            sdf_file.unlink()

@pytest.mark.parametrize("smiles", ["CCO", "O", "C1=CC=CC=C1"])
def test_sdf_file_creation(smiles, cleanup_sdf):
    directory = "tests/sdf_files"
    Path(directory).mkdir(exist_ok=True)

    filepath = sdf(smiles, directory, overwrite=True)

    assert Path(filepath).exists(), f"SDF file {filepath} was not created."

    supplier = Chem.SDMolSupplier(filepath)
    mol = supplier[0]
    assert mol is not None, "Failed to read molecule from SDF"
    assert mol.GetProp("SMILES") == smiles
