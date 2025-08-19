#!/usr/bin/env python3

import tempfile
from pathlib import Path
import unittest

import pytest

import screen
import util
import rcsb_pdb
import chem


# @pytest.fixture
def local_data_dir():
    return util.local_data_dir(__file__)


def test_dock_fabp4_with_bms309403_from_smiles(local_data_dir):
    receptor = util.mkdirs(screen.local_data_dir / "receptors") / "3RZY.pdb"
    output_dir = util.mkdirs(local_data_dir / "outputs")

    def dock_to_fabp4(compound, output_stem):
        return screen.dock(receptor, compound, output_dir / output_stem)

    smiles_docking = dock_to_fabp4(
        chem.sdf(
            "BMS",
            "[O-][S](=O)(=O)c1cccc2cccc(Nc3ccccc3)c12",
            directory_path=local_data_dir / "compounds",
        ),
        "bms_from_smiles_with_3rzy",
    )
    rcsb_docking = dock_to_fabp4(
        screen.pdb_partition(rcsb_pdb.write_pdb("2NNQ"), "T4B"),
        "bms_from_rcsb_with_3rzy",
    )
    assert util.files_are_equivalent(smiles_docking["sdf"], rcsb_docking["sdf"])


if __name__ == "__main__":
    test_dock_fabp4_with_bms309403_from_smiles(local_data_dir())
