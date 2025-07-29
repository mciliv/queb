import os
import argparse
import sys
import logging
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import SDWriter, AllChem

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)

def sdf(smiles: str, directory_path: str = ".", overwrite: bool = False) -> str:
    """Generates an SDF file from a SMILES string."""
    destination = Path(directory_path) / f"{smiles}.sdf"

    if not overwrite and destination.exists():
        logging.info(f"SDF already exists: {destination}")
        return str(destination)

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logging.error(f"Invalid SMILES string: {smiles}")
            return ""
            
        mol.SetProp("SMILES", smiles)
        Chem.Kekulize(mol)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)

        with SDWriter(str(destination)) as writer:
            writer.write(mol)

        logging.info(f"SDF file saved: {destination}")
        return str(destination)

    except Exception as e:
        logging.error(f"Error generating SDF for {smiles}: {e}")
        debug()
        return ""

def debug():
    """Debugging function to attach a debugger."""
    if os.getenv("PY_DEBUG") == "1":
        try:
            import debugpy
            print("Waiting for debugger to attach...")
            debugpy.listen(5678)
            debugpy.wait_for_client()
        except ImportError:
            print("debugpy not installed, skipping debug attach")

def main():
    parser = argparse.ArgumentParser(description="Generate an SDF file from a SMILES string.")
    parser.add_argument("smiles", type=str, help="The SMILES string to convert.")
    parser.add_argument("--dir", type=str, default=".", help="The directory to save the SDF file.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite the SDF file if it already exists.")

    args = parser.parse_args()

    try:
        sdf(args.smiles, args.dir, args.overwrite)
    except Exception as e:
        print(f"Exception occurred: {e}")
        raise

if __name__ == "__main__":
    main()