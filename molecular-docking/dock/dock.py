#!/usr/bin/env python3

import sys
from typing import List
import subprocess
import shutil
import argparse
import logging
from pathlib import Path
import time
from functools import partial

from Bio.PDB import PDBList

import filing
import processing
import drive
import log
import chem
import rcsb_pdb


def main():
    log.setup(Path.cwd())
    args = DockerArgumentParser().parse_args()
    Docker(args).dock_batch()

class DockerArgumentParser:
    def __init__(self):
        self.parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        self.parser.add_argument("-l", "--ligands", default=set(), help="Specify absolute path or a name", nargs="*", type=str)
        self.parser.add_argument("-f", "--default-ligands", action='store_true')
        self.parser.add_argument("-w", "--write-ligands", help="Overwrites existing ligands in the data_dir", action="store_true")
        self.parser.add_argument(
                "-e", "--ligand-excludes", help="Ligand file stems to exclude",
                default={"Insulin"}, nargs="*", type=str)
        self.parser.add_argument(
                "-b", "--complexes", default = {"2nnq", "5hz5"},
                nargs="*", type=str,
                help="RCSB PDB ids to extract receptor and use  for autobox")
        self.parser.add_argument("-p", "--apoproteins", default=set(), nargs="*", type=str, help="Gene names to search and use matching entries from RCSB PDB database")
        self.parser.add_argument("-d", "--data-dir-path", default=filing.local_data_dir(__file__))
        self.parser.add_argument("-o", "--dock-dir-name", default="docks")
        self.parser.add_argument("-s", "--sec-limit", default=float('inf'), help="Time limit for each dock", type=float)
        self.parser.add_argument("-a", "--redos", default=set(), nargs="*", type=str, help="List of output sdf file names to redo?")
        self.parser.add_argument("-g", "--gnina-options", nargs=argparse.REMAINDER, default=[], type=str, help="Override or include additional gnina option")
        self.args = self.parser.parse_args(args=[])
    
    def update_args(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self.args, key, value)

    def parse_args(self, args=sys.argv[1:]):
        self.args = self.parser.parse_args(args, namespace=self.args)
        return self.args

class Paths:
    def __init__(self, args):
        for attr in ["complexes", "receptors", "ligands", "docks"]:
            if attr == "docks":
                dir_name = args.dock_dir_name
            else: 
                dir_name = attr
            dir_path = args.data_dir_path / dir_name
            if not dir_path.exists():
                filing.renew(dir_name, mkdirs=True)
            setattr(self, attr, dir_path)

class Docker:
    def __init__(self, args: DockerArgumentParser):
        self.args = args
        self.paths = Paths(args)

        self.receptor_paths = set()
        self.ligand_paths = set()
        self._from_complexes()
        self._setup_apoproteins()
        self._setup_ligands()

    def _from_complexes(self):
        for complex_rcsb_pdb_id in self.args.complexes:
            complex_path = rcsb_pdb.write_pdb(complex_rcsb_pdb_id.lower(), self.paths.complexes)
            self.receptor_paths.add(chem.extract_from_structure(complex_path, chem.ProteinSelect(), self.paths.receptors / (complex_rcsb_pdb_id + chem.COMPLEX_PROTEIN_FILE_ENDING)))
            self.ligand_paths.add(chem.extract_from_structure(complex_path, chem.LigandSelect('A', rcsb_pdb.REF_LIGANDS[complex_rcsb_pdb_id]), self.paths.ligands / (complex_rcsb_pdb_id + chem.COMPLEX_LIGAND_FILE_ENDING)))

    def _setup_apoproteins(self):
         for apoprotein in self.args.apoproteins:
            rcsb_pdb.write_pdb(search_protein(apoprotein), path=self.path.receptors)

    def _setup_ligands(self):
        if self.args.write_ligands:
            filing.renew(self.paths.ligands, mkdirs=True, clean=True)
            chem.compounds(directory_path=self.paths.ligands)
        ligands = set(self.args.ligands)
        if self.args.default_ligands:
            ligands |= set(self.paths.ligands.iterdir())
        ligands |= {Path(path).stem.split('_')[0] for path in self.args.redos}
        ligands = {ligand for ligand in ligands if ligand.stem not in self.args.ligand_excludes}
        ligands = {Path(str(ligand).replace(' ', '_')) for ligand in ligands}
        self.ligand_paths |= filing.prepend_path(self.paths.ligands, ligands)

    def dock_batch(self):
        for receptor_path in self.receptor_paths:
            receptors_dock_dir_path = filing.renew(self.paths.docks / receptor_path.stem, mkdirs=True)
            for ligand_path in self.ligand_paths:
                self.dock(receptor_path, ligand_path, receptors_dock_dir_path)
        return self.paths.docks

    def dock(
            self, receptor_path: Path, ligand_path: Path,
            destination_dir: Path=Path.cwd()):
        dock_result = {ext: Path(destination_dir / (receptor_path.stem + "_" + ligand_path.stem)).with_suffix("." + ext) for ext in ("sdf", "txt")}
        if not (dock_result["sdf"].is_file() and dock_result["txt"].is_file() and chem.is_valid_sdf_with_molecule(dock_result["sdf"])):
            try:
                for path in dock_result.values(): filing.renew(path, mkdirs=False, clean=True)
                with open(dock_result["txt"], "a") as dock_txt_file:

                    logging.info(f"Starting {receptor_path} and {ligand_path}")
                    docking = subprocess.Popen(self.gnina_configured(receptor_path, ligand_path, dock_result["sdf"]), stdout=dock_txt_file, stderr=dock_txt_file)
                    processing.limit(docking, self.args.sec_limit, dock_result["txt"])
                chem.transfer_smiles_attribute(ligand_path, dock_result["sdf"])
                return dock_result["sdf"]
            except subprocess.CalledProcessError as e:
                logging.exception(f"Error processing {receptor_path} and {ligand_path}: {e}")
            except Exception as e:
                logging.exception(f"Error docking {receptor_path} and {ligand_path}: {e}")
        return None


    def gnina_configured(self, receptor_path, ligand_path, out_path) -> List[str]:
        return self.gnina_command(
            shutil.which("gnina") or Path.home() / "cure/gnina", receptor_path, ligand_path, self.paths.ligands / rcsb_pdb.ref_ligand(receptor_path), out_path)


    def gnina_command(self, gnina_path: Path, receptor_path: Path, ligand_path: Path, autobox_ligand: Path, out_path: Path) -> List[str]:
        return [
            gnina_path,
            "-r", receptor_path,
            "-l", ligand_path,
            "--autobox_ligand", autobox_ligand,
            "--autobox_extend", "1",
            "-o", out_path,
            "--seed", "0",
            "--exhaustiveness", "64",
        ] + self.args.gnina_options


if __name__ == "__main__":
    main()

