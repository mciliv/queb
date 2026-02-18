from pathlib import Path
import logging
from typing import Dict, List

import numpy as np
import pandas as pd
from Bio.PDB import PDBIO, PDBParser, Select, Structure, Model, Chain, Atom
from Bio.PDB.Residue import Residue

import spreadsheet
import filing


COMPLEX_PROTEIN_FILE_ENDING = "_protein.pdb"
COMPLEX_LIGAND_FILE_ENDING = "_ligand.pdb"


def local_data_dir():
    return filing.local_data_dir(__file__)


def compounds_dataframe() -> pd.DataFrame:
    smiless_with_names_range = "A:B"
    names_smiless_value_range = spreadsheet.get_values(
        "1VRl8dghA5t1JiLEa47OWW-hr-Qbvn9yksB4Hrga19Ck", smiless_with_names_range
    )
    compounds_spreadsheet = names_smiless_value_range["values"]
    compounds_dataframe = pd.DataFrame(compounds_spreadsheet[4:], columns=compounds_spreadsheet[3])
    compounds_dataframe["Name"] = compounds_dataframe["Name"].replace('', np.nan) 
    compounds_dataframe["Name"] = compounds_dataframe["Name"].fillna("index_" + compounds_dataframe.index.to_series().astype(str))
    return compounds_dataframe


def compound_dir():
    return Path(__file__).parent / "data/compounds"


def compounds(directory_path=compound_dir(), overwrite=False):
    for _, name, smiles in compounds_dataframe().itertuples():
        sdf(name.replace(' ', '_'), smiles, overwrite, directory_path)
    return directory_path


class ProteinSelect(Select):
    def accept_residue(self, residue):
        return residue.get_full_id()[3][0] == ' '


class LigandSelect(Select):
    def __init__(self, chain_id, residue_name):
        self.chain_id = chain_id
        self.residue_name = residue_name

    def accept_chain(self, chain):
        return chain.get_id() == self.chain_id

    def accept_residue(self, residue):
        return residue.get_resname() == self.residue_name


def extract_from_structure(pdb_path: Path, select: Select, output_path: Path):
    structure = PDBParser().get_structure(pdb_path.stem, str(pdb_path))
    return write_structure(structure, output_path, select)


def get_pdb_id_from_file(pdb_path):
    with open(pdb_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('HEADER'):
                pdb_id = line.split()[3]
                return pdb_id


def extract_residue(pdb_path: Path, identifier: str, output_path: Path) -> Dict[str, List[Residue]]:
    structure = PDBParser().get_structure(pdb_path.stem, str(pdb_path))
    for residue in structure.get_residues():
        if residue.get_resname() == identifier.upper():
            write_residue(residue, output_path)
            return residue


def write_residue(residue: Residue, output_path: Path) -> Path:
    write_structure(add_recursively(*ligand_scaffold(), residue), output_path)
    return output_path


def write_structure(structure: Structure, output_path: Path, select: Select=None) -> Path:
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    pdb_io.save(str(output_path), select)
    return output_path


def ligand_scaffold():
    ligand_structure = Structure.Structure("Ligand Structure")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    return (ligand_structure, model, chain)


def add_recursively(structure, model, chain, residue, atom):
    if atom is not None:
        residue.add(atom)
    chain.add(residue)
    model.add(chain)
    structure.add(model)
    return structure


def accumulate_parts(parts, super_part=Residue((" ", 0, " "), "LIG", " ")):
    for part in parts:
        super_part.add(part)
    return super_part


def pdb_partition(pdb_path: Path, identifier: str):
    pdb_partition_path = (pdb_path.parent / identifier).with_suffix(".pdb")
    with open(pdb_partition_path, "w") as pdb_path_file:
        subprocess.run(["grep", identifier, pdb_path], stdout=pdb_path_file)
    return pdb_partition_path


def transfer_smiles_attribute(source_path: Path, destination_path: Path):
    source_supplier = Chem.SDMolSupplier(str(source_path))
    smiless = []
    for mol in source_supplier:
        if mol:
            smiless.append(mol.GetProp('SMILES'))
    first_non_none = next((item for item in smiless if item), None)
    if first_non_none:
        destination_mols = mols_with_none(destination_path)
        set_smiles_for_mols(destination_path, destination_mols, [first_non_none] * len(destination_mols))
    else:
        add_smiless(destination_path)


def mols_with_none(file_path: Path) -> List[Chem.SDMolSupplier]:
    try:
        return [mol for mol in Chem.SDMolSupplier(str(file_path))]
    except OSError as e:
        return []


def add_smiless(file_path: Path) -> Path:
    mols = mols_with_none(file_path)
    smiless = [Chem.MolToSmiles(mol) if mol else None for mol in mols]
    return set_smiles_for_mols(file_path, mols, smiless) 


def set_smiles_for_mols(file_path: Path, mols: List[Chem.SDMolSupplier], smiless: List[str]) -> Path:
    assert len(mols) == len(smiless)
    with Chem.SDWriter(str(file_path)) as writer:
        for i in range(len(mols)):
            mols[i].SetProp('SMILES', smiless[i])
            writer.write(mols[i])
    return file_path

