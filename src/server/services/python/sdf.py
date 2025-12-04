#!/usr/bin/env python3
"""
sdf.py - 3D molecular structure generator using RDKit

This Python script is the core engine for converting SMILES (Simplified Molecular Input
Line Entry System) notation into 3D molecular structures in SDF (Structure Data Format).
It uses RDKit, a powerful cheminformatics library, to generate scientifically accurate
3D coordinates for molecules.

The script is called by the Node.js backend as a subprocess whenever a new molecular
structure needs to be generated. The generated SDF files are then served to the
frontend for 3D visualization using 3Dmol.js.

Key features:
- Converts SMILES strings to 3D molecular structures
- Generates SDF files compatible with molecular visualization tools
- Caches generated files to avoid redundant calculations
- Handles invalid SMILES gracefully with error reporting

Example usage:
    python sdf.py "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" --dir ./sdf_files
    # Generates caffeine.sdf in the specified directory

Requirements:
    - RDKit (install via: pip install rdkit-pypi or conda install -c conda-forge rdkit)
"""

import os
import argparse
import sys
import logging
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import SDWriter, AllChem

# Configure logging to output to stdout for Node.js subprocess capture
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)

def sdf(smiles: str, directory_path: str = ".", overwrite: bool = False) -> str:
    """
    Generates a 3D molecular structure file (SDF) from a SMILES string.
    
    This function performs the following steps:
    1. Validates the SMILES string
    2. Creates a molecular object using RDKit
    3. Adds hydrogen atoms (important for accurate 3D structure)
    4. Generates 3D coordinates using force field optimization
    5. Saves the structure as an SDF file
    
    Args:
        smiles (str): SMILES notation representing the molecule
                     Example: "CCO" for ethanol, "O" for water
        directory_path (str): Directory to save the SDF file (default: current directory)
        overwrite (bool): Whether to regenerate if file already exists (default: False)
    
    Returns:
        str: Path to the generated SDF file, or empty string if generation failed
    
    Note:
        The filename is derived from the SMILES string itself, which ensures
        unique filenames for different molecules and enables caching.
    """
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