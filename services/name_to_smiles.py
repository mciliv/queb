#!/usr/bin/env python3
import argparse
import sys


def convert_with_openeye(name: str) -> str:
    try:
        from openeye import oechem  # type: ignore
    except Exception:
        return "NONE"

    try:
        mol = oechem.OEGraphMol()
        if oechem.OEParseIUPACName(mol, name):
            return oechem.OEMolToSmiles(mol) or "NONE"
        return "NONE"
    except Exception:
        return "NONE"


def convert_with_rdkit(name: str) -> str:
    try:
        from rdkit import Chem  # type: ignore
    except Exception:
        return "NONE"

    try:
        mol = getattr(Chem, "MolFromName", None)
        if mol is None:
            # Fallback: try MolFromSmiles on the input if RDKit lacks MolFromName
            return "NONE"
        m = mol(name)
        if m:
            return Chem.MolToSmiles(m, isomericSmiles=True) or "NONE"
        return "NONE"
    except Exception:
        return "NONE"


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert chemical name to SMILES")
    parser.add_argument("--toolkit", choices=["openeye", "rdkit"], required=True)
    parser.add_argument("--name", required=True)
    args = parser.parse_args()

    name = args.name
    if not isinstance(name, str) or len(name.strip()) == 0:
        print("NONE")
        return 0

    if args.toolkit == "openeye":
        out = convert_with_openeye(name)
    else:
        out = convert_with_rdkit(name)

    # Always print a single line result for the Node caller to parse
    try:
        sys.stdout.write((out or "NONE").strip() + "\n")
        sys.stdout.flush()
    except Exception:
        # As a last resort
        print("NONE")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())



