import os
import argparse
import sys
import logging
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import SDWriter

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)

# Common mineral crystal structures and their unit cell parameters
CRYSTAL_STRUCTURES = {
    # Carbonates
    'CaCO3': {
        'name': 'Calcite (Calcium Carbonate)',
        'space_group': 'R-3c',
        'lattice': 'hexagonal',
        'a': 4.99, 'b': 4.99, 'c': 17.06,
        'alpha': 90, 'beta': 90, 'gamma': 120,
        'atoms': [
            ('Ca', [0.0, 0.0, 0.0]),
            ('C', [0.0, 0.0, 0.25]),
            ('O', [0.257, 0.0, 0.25]),
            ('O', [-0.257, 0.0, 0.25]),
            ('O', [0.0, 0.257, 0.25])
        ]
    },
    
    # Silicates  
    'SiO2': {
        'name': 'Quartz (Silicon Dioxide)',
        'space_group': 'P3121',
        'lattice': 'hexagonal',
        'a': 4.91, 'b': 4.91, 'c': 5.40,
        'alpha': 90, 'beta': 90, 'gamma': 120,
        'atoms': [
            ('Si', [0.470, 0.0, 0.0]),
            ('O', [0.415, 0.272, 0.119]),
            ('O', [0.415, 0.272, -0.119])
        ]
    },
    
    # Oxides
    'Al2O3': {
        'name': 'Corundum (Aluminum Oxide)', 
        'space_group': 'R-3c',
        'lattice': 'hexagonal',
        'a': 4.76, 'b': 4.76, 'c': 12.99,
        'alpha': 90, 'beta': 90, 'gamma': 120,
        'atoms': [
            ('Al', [0.0, 0.0, 0.352]),
            ('Al', [0.0, 0.0, -0.352]),
            ('O', [0.306, 0.0, 0.25]),
            ('O', [-0.306, 0.0, 0.25]),
            ('O', [0.0, 0.306, 0.25])
        ]
    },
    
    # Sulfides
    'FeS2': {
        'name': 'Pyrite (Iron Disulfide)',
        'space_group': 'Pa-3',
        'lattice': 'cubic',
        'a': 5.42, 'b': 5.42, 'c': 5.42,
        'alpha': 90, 'beta': 90, 'gamma': 90,
        'atoms': [
            ('Fe', [0.0, 0.0, 0.0]),
            ('S', [0.385, 0.385, 0.385]),
            ('S', [-0.385, -0.385, 0.385])
        ]
    },
    
    # Halides
    'NaCl': {
        'name': 'Halite (Sodium Chloride)',
        'space_group': 'Fm-3m', 
        'lattice': 'cubic',
        'a': 5.64, 'b': 5.64, 'c': 5.64,
        'alpha': 90, 'beta': 90, 'gamma': 90,
        'atoms': [
            ('Na', [0.0, 0.0, 0.0]),
            ('Cl', [0.5, 0.0, 0.0]),
            ('Na', [0.5, 0.5, 0.0]),
            ('Cl', [0.0, 0.5, 0.0])
        ]
    },
    
    # Water (Ice Ih structure)
    'H2O': {
        'name': 'Water (Ice Ih)',
        'space_group': 'P63/mmc',
        'lattice': 'hexagonal',
        'a': 4.52, 'b': 4.52, 'c': 7.36,
        'alpha': 90, 'beta': 90, 'gamma': 120,
        'atoms': [
            ('O', [0.333, 0.667, 0.25]),
            ('O', [0.667, 0.333, 0.75]),
            ('H', [0.333, 0.667, 0.125]),
            ('H', [0.667, 0.333, 0.875]),
            ('H', [0.333, 0.667, 0.375]),
            ('H', [0.667, 0.333, 0.625])
        ]
    },
    
    # Pure Elements (Face-Centered Cubic)
    'Au': {
        'name': 'Gold (Face-Centered Cubic)',
        'space_group': 'Fm-3m',
        'lattice': 'cubic',
        'a': 4.08, 'b': 4.08, 'c': 4.08,
        'alpha': 90, 'beta': 90, 'gamma': 90,
        'atoms': [
            ('Au', [0.0, 0.0, 0.0]),
            ('Au', [0.5, 0.5, 0.0]),
            ('Au', [0.5, 0.0, 0.5]),
            ('Au', [0.0, 0.5, 0.5])
        ]
    },
    'Cu': {
        'name': 'Copper (Face-Centered Cubic)',
        'space_group': 'Fm-3m',
        'lattice': 'cubic',
        'a': 3.61, 'b': 3.61, 'c': 3.61,
        'alpha': 90, 'beta': 90, 'gamma': 90,
        'atoms': [
            ('Cu', [0.0, 0.0, 0.0]),
            ('Cu', [0.5, 0.5, 0.0]),
            ('Cu', [0.5, 0.0, 0.5]),
            ('Cu', [0.0, 0.5, 0.5])
        ]
    },
    'Fe': {
        'name': 'Iron (Body-Centered Cubic)',
        'space_group': 'Im-3m',
        'lattice': 'cubic',
        'a': 2.87, 'b': 2.87, 'c': 2.87,
        'alpha': 90, 'beta': 90, 'gamma': 90,
        'atoms': [
            ('Fe', [0.0, 0.0, 0.0]),
            ('Fe', [0.5, 0.5, 0.5])
        ]
    },
    'Al': {
        'name': 'Aluminum (Face-Centered Cubic)',
        'space_group': 'Fm-3m',
        'lattice': 'cubic',
        'a': 4.05, 'b': 4.05, 'c': 4.05,
        'alpha': 90, 'beta': 90, 'gamma': 90,
        'atoms': [
            ('Al', [0.0, 0.0, 0.0]),
            ('Al', [0.5, 0.5, 0.0]),
            ('Al', [0.5, 0.0, 0.5]),
            ('Al', [0.0, 0.5, 0.5])
        ]
    }
}

def parse_mineral_formula(formula):
    """Simple parser for common mineral formulas"""
    # Normalize formula (remove spaces, handle common variations)
    formula = formula.replace(' ', '').strip()
    
    # Direct lookup first
    if formula in CRYSTAL_STRUCTURES:
        return formula
        
    # Handle common variations
    variations = {
        'CaCO₃': 'CaCO3',
        'SiO₂': 'SiO2', 
        'Al₂O₃': 'Al2O3',
        'FeS₂': 'FeS2',
        'quartz': 'SiO2',
        'calcite': 'CaCO3',
        'corundum': 'Al2O3',
        'pyrite': 'FeS2',
        'halite': 'NaCl',
        'salt': 'NaCl'
    }
    
    return variations.get(formula.lower(), formula)

def generate_crystal_coordinates(crystal_data, repeat_x=2, repeat_y=2, repeat_z=2):
    """Generate 3D coordinates for crystal unit cell with repetitions"""
    coordinates = []
    
    # Unit cell parameters
    a, b, c = crystal_data['a'], crystal_data['b'], crystal_data['c'] 
    alpha, beta, gamma = np.radians([crystal_data['alpha'], crystal_data['beta'], crystal_data['gamma']])
    
    # Create transformation matrix for fractional to Cartesian coordinates
    cos_alpha, cos_beta, cos_gamma = np.cos(alpha), np.cos(beta), np.cos(gamma)
    sin_gamma = np.sin(gamma)
    
    # Volume calculation
    volume = a * b * c * np.sqrt(1 + 2*cos_alpha*cos_beta*cos_gamma - cos_alpha**2 - cos_beta**2 - cos_gamma**2)
    
    # Transformation matrix  
    transform = np.array([
        [a, b*cos_gamma, c*cos_beta],
        [0, b*sin_gamma, c*(cos_alpha - cos_beta*cos_gamma)/sin_gamma],
        [0, 0, volume/(a*b*sin_gamma)]
    ])
    
    atom_id = 1
    for x_rep in range(repeat_x):
        for y_rep in range(repeat_y):
            for z_rep in range(repeat_z):
                for element, frac_coords in crystal_data['atoms']:
                    # Apply unit cell translations
                    frac_x = frac_coords[0] + x_rep
                    frac_y = frac_coords[1] + y_rep  
                    frac_z = frac_coords[2] + z_rep
                    
                    # Convert to Cartesian coordinates
                    frac_vec = np.array([frac_x, frac_y, frac_z])
                    cart_coords = transform @ frac_vec
                    
                    coordinates.append({
                        'atom_id': atom_id,
                        'element': element,
                        'x': cart_coords[0],
                        'y': cart_coords[1], 
                        'z': cart_coords[2]
                    })
                    atom_id += 1
                    
    return coordinates

def create_crystal_sdf(formula, coordinates, crystal_data):
    """Create SDF content for crystal structure"""
    
    # SDF header
    sdf_content = f"""  {crystal_data['name']}
  Generated crystal structure
  
{len(coordinates):3d}  0  0  0  0  0  0  0  0  0999 V2000
"""
    
    # Atom block
    for coord in coordinates:
        sdf_content += f"{coord['x']:10.4f}{coord['y']:10.4f}{coord['z']:10.4f} {coord['element']:<3s} 0  0  0  0  0  0  0  0  0  0  0  0\n"
    
    # Properties block
    sdf_content += """M  END
> <FORMULA>
""" + formula + """

> <CRYSTAL_SYSTEM>
""" + crystal_data['lattice'] + """

> <SPACE_GROUP>
""" + crystal_data['space_group'] + """

> <MINERAL_NAME>
""" + crystal_data['name'] + """

$$$$
"""
    
    return sdf_content

def crystal(formula: str, directory_path: str = ".", overwrite: bool = False) -> str:
    """Generates an SDF file from a mineral formula."""
    
    # Parse and normalize the formula
    normalized_formula = parse_mineral_formula(formula)
    
    # Use normalized formula for filename, but original for display
    safe_filename = normalized_formula.replace('/', '_').replace('(', '_').replace(')', '_')
    destination = Path(directory_path) / f"{safe_filename}.sdf"

    if not overwrite and destination.exists():
        logging.info(f"Crystal SDF already exists: {destination}")
        return str(destination)

    try:
        # Look up crystal structure
        if normalized_formula not in CRYSTAL_STRUCTURES:
            logging.error(f"Unknown mineral formula: {formula} (normalized: {normalized_formula})")
            logging.info(f"Known minerals: {list(CRYSTAL_STRUCTURES.keys())}")
            return ""
            
        crystal_data = CRYSTAL_STRUCTURES[normalized_formula]
        logging.info(f"Generating crystal structure for {crystal_data['name']}")
        
        # Generate coordinates for a 2x2x2 unit cell repetition
        coordinates = generate_crystal_coordinates(crystal_data, repeat_x=2, repeat_y=2, repeat_z=2)
        
        # Create SDF content
        sdf_content = create_crystal_sdf(normalized_formula, coordinates, crystal_data)
        
        # Write to file
        with open(destination, 'w') as f:
            f.write(sdf_content)

        logging.info(f"Crystal SDF file saved: {destination}")
        return str(destination)

    except Exception as e:
        logging.error(f"Error generating crystal structure for {formula}: {e}")
        return ""

def main():
    parser = argparse.ArgumentParser(description="Generate a crystal SDF file from a mineral formula.")
    parser.add_argument("formula", type=str, help="The mineral formula to convert (e.g., CaCO3, SiO2).")
    parser.add_argument("--dir", type=str, default=".", help="The directory to save the SDF file.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite the SDF file if it already exists.")

    args = parser.parse_args()

    try:
        result = crystal(args.formula, args.dir, args.overwrite)
        if result:
            print(f"✅ Crystal structure generated: {result}")
        else:
            print(f"❌ Failed to generate crystal structure for: {args.formula}")
            sys.exit(1)
    except Exception as e:
        print(f"Exception occurred: {e}")
        raise

if __name__ == "__main__":
    main() 