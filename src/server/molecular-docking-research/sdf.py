#!/usr/bin/env python3
"""
SDF Retriever - Unified SDF file fetching utility

This module provides a clean interface for retrieving SDF files from:
1. Local file system (existing SDF directories)
2. Project's internal API endpoints
3. Direct PubChem API calls

Usage:
    python sdf.py caffeine --dir ./output
    python sdf.py CCO
    python sdf.py "calcium carbonate"
"""

import os
import argparse
import sys
import logging
from pathlib import Path
import urllib.request
import urllib.parse
import json

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)

# Common molecular formula variations for lookup
FORMULA_VARIATIONS = {
    'CaCO₃': 'CaCO3',
    'SiO₂': 'SiO2', 
    'Al₂O₃': 'Al2O3',
    'FeS₂': 'FeS2',
    'quartz': 'SiO2',
    'calcite': 'CaCO3',
    'corundum': 'Al2O3',
    'pyrite': 'FeS2',
    'halite': 'NaCl',
    'salt': 'NaCl',
    'water': 'H2O',
    'ice': 'H2O',
    'ethanol': 'CCO',
    'alcohol': 'CCO',
    'caffeine': 'caffeine',
    'aspirin': 'aspirin',
    'glucose': 'glucose'
}

def normalize_formula(formula):
    """Normalize molecular formula for consistent lookup"""
    formula = formula.replace(' ', '').strip()
    return FORMULA_VARIATIONS.get(formula.lower(), formula)

def find_local_sdf(formula, project_root=None):
    """Search for existing SDF files in project directories"""
    if not project_root:
        project_root = Path(__file__).parent.parent.parent
    
    # Define SDF directories to search
    sdf_directories = [
        project_root / "backend" / "sdf_files",
        project_root / "tests" / "sdf_files"
    ]
    
    # Generate possible filename variations
    normalized = normalize_formula(formula)
    possible_filenames = [
        f"{normalized}.sdf",
        f"_{normalized.replace('2', '_2_').replace('3', '_3_')}.sdf",
        f"{normalized.lower()}.sdf",
        f"{formula.lower()}.sdf"
    ]
    
    # Search for exact matches first
    for sdf_dir in sdf_directories:
        if not sdf_dir.exists():
            continue
            
        for filename in possible_filenames:
            sdf_path = sdf_dir / filename
            if sdf_path.exists():
                logging.info(f"Found local SDF: {sdf_path}")
                return str(sdf_path)
    
    # Search for partial matches
    for sdf_dir in sdf_directories:
        if not sdf_dir.exists():
            continue
            
        try:
            for sdf_file in sdf_dir.glob("*.sdf"):
                filename_lower = sdf_file.stem.lower()
                formula_lower = normalized.lower()
                
                # Check for partial matches
                if (formula_lower in filename_lower or 
                    filename_lower in formula_lower or
                    any(part in filename_lower for part in formula_lower.split() if len(part) > 2)):
                    
                    logging.info(f"Found partial match: {sdf_file}")
                    return str(sdf_file)
        except Exception as e:
            logging.debug(f"Error searching {sdf_dir}: {e}")
    
    return None

def fetch_via_internal_api(formula, api_base="http://localhost:8080"):
    """Fetch SDF via the project's internal API"""
    try:
        api_url = f"{api_base}/pubchem/sdf"
        payload = json.dumps({"name": formula}).encode('utf-8')
        req = urllib.request.Request(api_url, data=payload, method='POST')
        req.add_header('Content-Type', 'application/json')

        with urllib.request.urlopen(req, timeout=15) as response:
            if response.status == 200:
                response_data = json.loads(response.read().decode('utf-8'))
                if response_data.get('sdf'):
                    logging.info("Retrieved SDF via internal API")
                    return response_data['sdf']
                logging.debug("API response missing SDF data")
            else:
                logging.debug(f"API returned status {response.status}")

    except Exception as e:
        logging.debug(f"Internal API call failed: {e}")

    return None

def fetch_via_pubchem(formula):
    """Fetch SDF directly from PubChem API"""
    try:
        # Try compound search by name
        encoded_name = urllib.parse.quote(formula)
        pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/SDF"
        
        with urllib.request.urlopen(pubchem_url, timeout=15) as response:
            if response.status == 200:
                sdf_content = response.read().decode('utf-8')
                if len(sdf_content) > 100:  # Basic validation
                    logging.info(f"Retrieved SDF from PubChem")
                    return sdf_content
    
    except Exception as e:
        logging.debug(f"PubChem API call failed: {e}")
    
    return None

def validate_sdf(content):
    """Basic SDF format validation"""
    if not content or len(content) < 50:
        return False
    
    # Check for SDF format markers
    has_end_marker = "M  END" in content or "$$$$" in content
    has_reasonable_length = len(content) > 100
    
    return has_end_marker and has_reasonable_length

def retrieve_sdf(formula, output_dir=".", overwrite=False):
    """
    Main SDF retrieval function with multiple fallback strategies
    
    Args:
        formula: Molecular formula or compound name
        output_dir: Directory to save the SDF file
        overwrite: Whether to overwrite existing files
    
    Returns:
        Path to the SDF file if successful, None otherwise
    """
    normalized = normalize_formula(formula)
    output_path = Path(output_dir) / f"{normalized}.sdf"
    
    # Check if output file already exists
    if not overwrite and output_path.exists():
        logging.info(f"SDF already exists: {output_path}")
        return str(output_path), False
    
    # Strategy 1: Look for local files
    local_sdf = find_local_sdf(formula)
    if local_sdf:
        if str(output_path) != local_sdf:
            # Copy to output directory
            try:
                output_path.parent.mkdir(parents=True, exist_ok=True)
                import shutil
                shutil.copy2(local_sdf, output_path)
                logging.info(f"Copied local SDF to: {output_path}")
                return str(output_path), True
            except Exception as e:
                logging.error(f"Failed to copy SDF: {e}")
                return local_sdf, False  # Return original path
        else:
            return local_sdf, False
    
    # Strategy 2: Try internal API (if server is running)
    sdf_content = fetch_via_internal_api(formula)
    if sdf_content and validate_sdf(sdf_content):
        try:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, 'w') as f:
                f.write(sdf_content)
            logging.info(f"Saved SDF from internal API: {output_path}")
            return str(output_path), True
        except Exception as e:
            logging.error(f"Failed to save SDF: {e}")
    
    # Strategy 3: Direct PubChem API
    sdf_content = fetch_via_pubchem(formula)
    if sdf_content and validate_sdf(sdf_content):
        try:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, 'w') as f:
                f.write(sdf_content)
            logging.info(f"Saved SDF from PubChem: {output_path}")
            return str(output_path), True
        except Exception as e:
            logging.error(f"Failed to save SDF: {e}")
    
    # All strategies failed
    logging.error(f"Could not retrieve SDF for: {formula}")
    return None, False

def get_solution_info(result_path):
    """Determine which solution strategy was used"""
    if not result_path:
        return "❌ No solution found"
    
    path_str = str(result_path)
    if "backend/sdf_files" in path_str or "tests/sdf_files" in path_str:
        return "📁 Local file system"
    else:
        return "🌐 API download"

def main():
    parser = argparse.ArgumentParser(
        description="Retrieve SDF files for molecular compounds",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python sdf.py caffeine
  python sdf.py CCO --dir ./molecules
  python sdf.py "calcium carbonate" --overwrite
        """
    )
    
    parser.add_argument("formula", 
                       help="Molecular formula or compound name (e.g., CCO, caffeine, CaCO3)")
    parser.add_argument("--dir", default=".", 
                       help="Directory to save the SDF file (default: current directory)")
    parser.add_argument("--overwrite", action="store_true",
                       help="Overwrite existing SDF files")
    parser.add_argument("--quiet", action="store_true",
                       help="Suppress info messages")
    parser.add_argument("--open", action="store_true",
                       help="Open the SDF file after retrieval using the system viewer")
    parser.add_argument("--app", default=None,
                       help="Specific app to open with (e.g., 'Avogadro2'). macOS only.")
    parser.add_argument("--publish", action="store_true",
                       help="Copy SDF into the app's served sdf_files directory for web viewing")
    parser.add_argument("--open-viewer", action="store_true",
                       help="Open the SDF in the app's web viewer (/sdf_files/<file>)")
    parser.add_argument("--base-url", default="http://localhost:8080",
                       help="Base URL of the running app for viewer open (default: http://localhost:8080)")
    
    args = parser.parse_args()
    
    if args.quiet:
        logging.getLogger().setLevel(logging.WARNING)
    
    try:
        result_path, changed = retrieve_sdf(args.formula, args.dir, args.overwrite)
        
        if result_path:
            print(f"✅ SDF file retrieved: {result_path}")
            print(f"   {get_solution_info(result_path)}")
            
            # Show file info
            try:
                file_size = Path(result_path).stat().st_size
                print(f"   📄 File size: {file_size:,} bytes")
            except:
                pass

            # Optionally publish into served directory for web viewer
            published_path = None
            published_url = None
            if args.publish or args.open_viewer:
                try:
                    project_root = Path(__file__).parent.parent.parent
                    node_env = os.environ.get('NODE_ENV', '').lower()
                    served_dir = project_root / ("tests/sdf_files" if node_env == 'test' else "backend/sdf_files")
                    served_dir.mkdir(parents=True, exist_ok=True)
                    published_path = served_dir / Path(result_path).name
                    if str(published_path) != str(result_path):
                        import shutil
                        shutil.copy2(result_path, published_path)
                    published_url = f"{args.base_url}/sdf_files/{published_path.name}"
                    print(f"   🌐 Published for web viewer: {published_url}")
                except Exception as e:
                    print(f"   ⚠️  Could not publish to served directory: {e}")

            # Optionally open the file (native app)
            if args.open:
                try:
                    import platform
                    import subprocess
                    system = platform.system().lower()
                    if system == 'darwin':
                        # macOS
                        if changed:
                            # New/updated file: always open
                            if args.app:
                                subprocess.run(["open", "-a", args.app, result_path], check=True)
                            else:
                                subprocess.run(["open", result_path], check=True)
                        else:
                            # Existing file: open only if explicit app specified
                            if args.app:
                                subprocess.run(["open", "-a", args.app, result_path], check=True)
                        print("   🚀 Opened in viewer")
                    elif system == 'windows':
                        if changed or args.app:
                            os.startfile(result_path)  # type: ignore[attr-defined]
                        print("   🚀 Opened in default viewer")
                    else:
                        # Linux and others
                        # Try xdg-open
                        if changed or args.app:
                            subprocess.run(["xdg-open", result_path], check=True)
                        print("   🚀 Opened in default viewer")
                except Exception as e:
                    print(f"   ⚠️  Could not open viewer automatically: {e}")
                    if platform.system().lower() == 'darwin':
                        print("   Tip: install Avogadro2: brew install --cask avogadro2")
                        print("   Then run: open -a Avogadro2 " + result_path)

            # Optionally open in our web viewer
            if args.open_viewer:
                try:
                    import webbrowser
                    target_url = published_url
                    if not target_url:
                        # Best-effort fallback to raw file name
                        target_url = f"{args.base_url}/sdf_files/{Path(result_path).name}"
                    webbrowser.open(target_url)
                    print("   🌐 Opened in app viewer: " + target_url)
                except Exception as e:
                    print(f"   ⚠️  Could not open app viewer: {e}")
        else:
            print(f"❌ No SDF file found for: {args.formula}")
            print(f"   Try checking the formula spelling or network connection")
            sys.exit(1)
    
    except KeyboardInterrupt:
        print("\n⚠️  Operation cancelled by user")
        sys.exit(1)
    except Exception as e:
        print(f"💥 Error: {e}")
        if not args.quiet:
            raise
        sys.exit(1)

if __name__ == "__main__":
    main()
