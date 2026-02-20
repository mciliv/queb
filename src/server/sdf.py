#!/usr/bin/env python3
"""
sdf.py - Unified SDF retrieval and generation

Pipeline (in order):
1. Local cache
2. PubChem (by name or SMILES)
3. RDKit 3D generation (requires rdkit)

Usage:
    python sdf.py CCO --dir ./sdf_files
    python sdf.py caffeine --dir ./sdf_files
    python sdf.py "calcium carbonate" --overwrite
"""

import os
import sys
import argparse
import logging
import shutil
import urllib.request
import urllib.parse
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)

FORMULA_VARIATIONS = {
    'water': 'H2O', 'ice': 'H2O',
    'ethanol': 'CCO', 'alcohol': 'CCO',
    'salt': 'NaCl', 'halite': 'NaCl',
    'quartz': 'SiO2', 'calcite': 'CaCO3',
    'corundum': 'Al2O3', 'pyrite': 'FeS2',
}

def normalize(formula):
    return FORMULA_VARIATIONS.get(formula.strip().lower(), formula.strip())

def safe_filename(formula):
    import hashlib
    cleaned = ''.join(c if c.isalnum() or c in '._-' else ('__' if c == '=' else '_') for c in formula)
    if len(cleaned) > 50:
        h = hashlib.md5(formula.encode()).hexdigest()[:12]
        return f"{cleaned[:8]}_{h}.sdf"
    return f"{cleaned}.sdf"

def find_local(formula, output_dir):
    name = normalize(formula)
    root = Path(__file__).parent.parent.parent
    dirs = [Path(output_dir), root / "public" / "sdf_files", root / "tests" / "sdf_files"]
    for d in dirs:
        for filename in [f"{name}.sdf", safe_filename(name)]:
            p = d / filename
            if p.exists() and p.stat().st_size > 50:
                return str(p)
    return None

def fetch_pubchem(formula):
    encoded = urllib.parse.quote(formula)
    for endpoint in [
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/SDF",
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded}/SDF",
    ]:
        try:
            with urllib.request.urlopen(endpoint, timeout=15) as r:
                if r.status == 200:
                    content = r.read().decode('utf-8')
                    if len(content) > 100 and ("M  END" in content or "$$$$" in content):
                        return content
        except Exception:
            continue
    return None

def generate_rdkit(smiles, dest):
    from rdkit import Chem
    from rdkit.Chem import SDWriter, AllChem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mol.SetProp("SMILES", smiles)
    Chem.Kekulize(mol)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    with SDWriter(str(dest)) as w:
        w.write(mol)
    return True

def sdf(formula, output_dir=".", overwrite=False):
    """
    Retrieve or generate an SDF file for a compound name or SMILES string.
    Returns path to the SDF file, or empty string on failure.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    name = normalize(formula)
    dest = output_dir / safe_filename(name)

    # 1. Local cache
    if not overwrite:
        local = find_local(formula, output_dir)
        if local:
            logging.info(f"SDF cached: {local}")
            return local

    # 2. PubChem
    content = fetch_pubchem(name)
    if content:
        dest.write_text(content)
        logging.info(f"SDF from PubChem: {dest}")
        return str(dest)

    # 3. RDKit (SMILES input)
    try:
        if generate_rdkit(name, dest):
            logging.info(f"SDF from RDKit: {dest}")
            return str(dest)
    except ImportError:
        logging.debug("rdkit not available")
    except Exception as e:
        logging.error(f"RDKit failed: {e}")

    logging.error(f"Could not retrieve SDF for: {formula}")
    return ""

def debug():
    if os.getenv("PY_DEBUG") == "1":
        try:
            import debugpy
            print("Waiting for debugger...")
            debugpy.listen(5678)
            debugpy.wait_for_client()
        except ImportError:
            pass

def main():
    parser = argparse.ArgumentParser(
        description="Retrieve or generate SDF files for molecular compounds"
    )
    parser.add_argument("formula", help="SMILES string or compound name (e.g., CCO, caffeine)")
    parser.add_argument("--dir", default=".", help="Output directory (default: .)")
    parser.add_argument("--overwrite", action="store_true", help="Regenerate even if cached")
    parser.add_argument("--quiet", action="store_true", help="Suppress info messages")
    parser.add_argument("--open", action="store_true", help="Open in system viewer after retrieval")
    parser.add_argument("--app", default=None, help="App to open with (macOS, e.g. 'Avogadro2')")
    parser.add_argument("--publish", action="store_true", help="Copy into served sdf_files dir")
    parser.add_argument("--open-viewer", action="store_true", help="Open in app web viewer")
    parser.add_argument("--base-url", default="http://localhost:8080", help="App base URL for viewer")
    args = parser.parse_args()

    if args.quiet:
        logging.getLogger().setLevel(logging.WARNING)

    try:
        result = sdf(args.formula, args.dir, args.overwrite)

        if not result:
            print(f"❌ No SDF file found for: {args.formula}")
            sys.exit(1)

        result_path = Path(result)
        print(f"✅ SDF file retrieved: {result}")
        print(f"   📄 File size: {result_path.stat().st_size:,} bytes")

        published_url = None
        if args.publish or args.open_viewer:
            try:
                root = Path(__file__).parent.parent.parent
                served = root / ("tests/sdf_files" if os.environ.get("NODE_ENV") == "test" else "public/sdf_files")
                served.mkdir(parents=True, exist_ok=True)
                published = served / result_path.name
                if str(published) != str(result_path):
                    shutil.copy2(result, published)
                published_url = f"{args.base_url}/sdf_files/{result_path.name}"
                print(f"   🌐 Published: {published_url}")
            except Exception as e:
                print(f"   ⚠️  Could not publish: {e}")

        if args.open:
            try:
                import platform, subprocess
                if platform.system() == 'Darwin':
                    cmd = ["open", "-a", args.app, result] if args.app else ["open", result]
                    subprocess.run(cmd, check=True)
                elif platform.system() == 'Windows':
                    os.startfile(result)
                else:
                    subprocess.run(["xdg-open", result], check=True)
                print("   🚀 Opened in viewer")
            except Exception as e:
                print(f"   ⚠️  Could not open viewer: {e}")

        if args.open_viewer:
            try:
                import webbrowser
                url = published_url or f"{args.base_url}/sdf_files/{result_path.name}"
                webbrowser.open(url)
                print(f"   🌐 Opened in web viewer: {url}")
            except Exception as e:
                print(f"   ⚠️  Could not open web viewer: {e}")

    except KeyboardInterrupt:
        print("\n⚠️  Cancelled")
        sys.exit(1)
    except Exception as e:
        debug()
        print(f"💥 Error: {e}")
        if not args.quiet:
            raise
        sys.exit(1)

if __name__ == "__main__":
    main()
