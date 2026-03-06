const fs = require('fs');
const path = require('path');

// Van der Waals radii in Angstroms (for rendering)
const VDW_RADII = {
  H: 1.20, He: 1.40, Li: 1.82, Be: 1.53, B: 1.92, C: 1.70, N: 1.55, O: 1.52,
  F: 1.47, Ne: 1.54, Na: 2.27, Mg: 1.73, Al: 1.84, Si: 2.10, P: 1.80, S: 1.80,
  Cl: 1.75, Ar: 1.88, K: 2.75, Ca: 2.31, Fe: 2.04, Cu: 1.40, Zn: 1.39,
  Br: 1.85, I: 1.98, Au: 1.66, Ag: 1.72
};

// CPK colors (hex)
const CPK_COLORS = {
  H: '#FFFFFF', He: '#D9FFFF', Li: '#CC80FF', Be: '#C2FF00', B: '#FFB5B5',
  C: '#909090', N: '#3050F8', O: '#FF0D0D', F: '#90E050', Ne: '#B3E3F5',
  Na: '#AB5CF2', Mg: '#8AFF00', Al: '#BFA6A6', Si: '#F0C8A0', P: '#FF8000',
  S: '#FFFF30', Cl: '#1FF01F', Ar: '#80D1E3', K: '#8F40D4', Ca: '#3DFF00',
  Fe: '#E06633', Cu: '#C88033', Zn: '#7D80B0', Br: '#A62929', I: '#940094',
  Au: '#FFD123', Ag: '#C0C0C0'
};

class MaterialSceneBuilder {
  constructor(dependencies = {}) {
    this.logger = dependencies.logger || console;
    this._fixtureDir = path.join(__dirname, '..', '..', '..', 'resources', 'material-scenes');
  }

  build(materialData) {
    const mode = materialData.visualization_mode;
    switch (mode) {
      case 'crystal':
        return this._buildCrystal(materialData);
      case 'liquid':
        return this._buildLiquid(materialData);
      case 'gas':
        return this._buildGas(materialData);
      default:
        return null;
    }
  }

  _buildCrystal(materialData) {
    const md = materialData.material_data || {};
    const latticeType = md.lattice_type;
    const lc = md.lattice_constants || { a: 4.0, b: 4.0, c: 4.0 };

    let unitCellAtoms;
    if (md.basis_atoms && md.basis_atoms.length > 0) {
      unitCellAtoms = md.basis_atoms;
    } else {
      unitCellAtoms = this._getDefaultBasis(latticeType, materialData);
    }

    // Replicate unit cell in a 3x3x3 supercell
    const repeats = 3;
    const atoms = [];
    for (let ix = 0; ix < repeats; ix++) {
      for (let iy = 0; iy < repeats; iy++) {
        for (let iz = 0; iz < repeats; iz++) {
          for (const basis of unitCellAtoms) {
            const x = (basis.x + ix) * lc.a;
            const y = (basis.y + iy) * (lc.b || lc.a);
            const z = (basis.z + iz) * (lc.c || lc.a);
            atoms.push({
              element: basis.element,
              x, y, z,
              radius: VDW_RADII[basis.element] || 1.5,
              color: CPK_COLORS[basis.element] || '#FF69B4'
            });
          }
        }
      }
    }

    // Center the supercell around origin
    const cx = (repeats * lc.a) / 2;
    const cy = (repeats * (lc.b || lc.a)) / 2;
    const cz = (repeats * (lc.c || lc.a)) / 2;
    for (const atom of atoms) {
      atom.x -= cx;
      atom.y -= cy;
      atom.z -= cz;
    }

    return {
      atoms,
      box: {
        a: repeats * lc.a,
        b: repeats * (lc.b || lc.a),
        c: repeats * (lc.c || lc.a)
      },
      metadata: {
        visualization_mode: 'crystal',
        lattice_type: latticeType,
        space_group: md.space_group,
        description: materialData.description,
        atom_count: atoms.length
      }
    };
  }

  _getDefaultBasis(latticeType, materialData) {
    switch (latticeType) {
      case 'FCC':
        return this._fccBasis(materialData);
      case 'BCC':
        return this._bccBasis(materialData);
      case 'diamond_cubic':
        return this._diamondCubicBasis(materialData);
      case 'rock_salt':
        return this._rockSaltBasis(materialData);
      case 'HCP':
        return this._hcpBasis(materialData);
      case 'hexagonal_ice':
        return this._hexagonalIceBasis();
      case 'simple_cubic':
        return this._simpleCubicBasis(materialData);
      default:
        return this._simpleCubicBasis(materialData);
    }
  }

  _simpleCubicBasis(materialData) {
    const el = this._primaryElement(materialData);
    return [{ element: el, x: 0, y: 0, z: 0 }];
  }

  _bccBasis(materialData) {
    const el = this._primaryElement(materialData);
    return [
      { element: el, x: 0, y: 0, z: 0 },
      { element: el, x: 0.5, y: 0.5, z: 0.5 }
    ];
  }

  _fccBasis(materialData) {
    const el = this._primaryElement(materialData);
    return [
      { element: el, x: 0, y: 0, z: 0 },
      { element: el, x: 0.5, y: 0.5, z: 0 },
      { element: el, x: 0.5, y: 0, z: 0.5 },
      { element: el, x: 0, y: 0.5, z: 0.5 }
    ];
  }

  _diamondCubicBasis(materialData) {
    const el = this._primaryElement(materialData);
    return [
      { element: el, x: 0, y: 0, z: 0 },
      { element: el, x: 0.5, y: 0.5, z: 0 },
      { element: el, x: 0.5, y: 0, z: 0.5 },
      { element: el, x: 0, y: 0.5, z: 0.5 },
      { element: el, x: 0.25, y: 0.25, z: 0.25 },
      { element: el, x: 0.75, y: 0.75, z: 0.25 },
      { element: el, x: 0.75, y: 0.25, z: 0.75 },
      { element: el, x: 0.25, y: 0.75, z: 0.75 }
    ];
  }

  _rockSaltBasis(materialData) {
    const chemicals = materialData.chemicals || [];
    let cation = 'Na', anion = 'Cl';
    if (chemicals.length > 0) {
      const name = (chemicals[0].name || '').toLowerCase();
      if (name.includes('sodium') || name.includes('nacl')) {
        cation = 'Na'; anion = 'Cl';
      } else if (name.includes('potassium') || name.includes('kcl')) {
        cation = 'K'; anion = 'Cl';
      } else if (name.includes('lithium') || name.includes('lif')) {
        cation = 'Li'; anion = 'F';
      }
    }
    return [
      { element: cation, x: 0, y: 0, z: 0 },
      { element: anion, x: 0.5, y: 0, z: 0 },
      { element: cation, x: 0.5, y: 0.5, z: 0 },
      { element: anion, x: 0, y: 0.5, z: 0 },
      { element: cation, x: 0, y: 0, z: 0.5 },
      { element: anion, x: 0.5, y: 0, z: 0.5 },
      { element: cation, x: 0.5, y: 0.5, z: 0.5 },
      { element: anion, x: 0, y: 0.5, z: 0.5 }
    ];
  }

  _hcpBasis(materialData) {
    const el = this._primaryElement(materialData);
    return [
      { element: el, x: 0, y: 0, z: 0 },
      { element: el, x: 1/3, y: 2/3, z: 0.5 }
    ];
  }

  _hexagonalIceBasis() {
    // Simplified ice Ih: oxygen positions in wurtzite-like arrangement
    // with hydrogens placed along O-O bonds
    return [
      { element: 'O', x: 0, y: 0, z: 0 },
      { element: 'H', x: 0.1, y: 0, z: 0.06 },
      { element: 'H', x: -0.05, y: 0.09, z: 0.06 },
      { element: 'O', x: 1/3, y: 2/3, z: 0.5 },
      { element: 'H', x: 0.43, y: 0.67, z: 0.56 },
      { element: 'H', x: 0.28, y: 0.76, z: 0.56 }
    ];
  }

  _primaryElement(materialData) {
    const chemicals = materialData.chemicals || [];
    if (chemicals.length > 0) {
      const name = (chemicals[0].name || '').toLowerCase();
      if (name.includes('iron') || name === 'fe') return 'Fe';
      if (name.includes('gold') || name === 'au') return 'Au';
      if (name.includes('silver') || name === 'ag') return 'Ag';
      if (name.includes('copper') || name === 'cu') return 'Cu';
      if (name.includes('aluminum') || name.includes('aluminium')) return 'Al';
      if (name.includes('carbon') || name.includes('diamond')) return 'C';
      if (name.includes('silicon')) return 'Si';
    }
    return 'Fe';
  }

  _buildLiquid(materialData) {
    // Try to load pre-computed fixture
    const md = materialData.material_data || {};
    const formula = (md.molecular_formula || '').toLowerCase();

    if (formula === 'h2o' || (materialData.object || '').toLowerCase().includes('water')) {
      return this._loadFixture('water-liquid-298K.json');
    }

    // For unknown liquids, generate a simple random box
    return this._generateRandomLiquidBox(materialData);
  }

  _loadFixture(filename) {
    const filePath = path.join(this._fixtureDir, filename);
    try {
      const data = JSON.parse(fs.readFileSync(filePath, 'utf8'));
      return data;
    } catch (err) {
      this.logger.warn && this.logger.warn(`Failed to load fixture ${filename}: ${err.message}`);
      return null;
    }
  }

  _generateRandomLiquidBox(materialData) {
    // Fallback: generate a random arrangement
    // This is a placeholder - real equilibrated configs would come from MD simulations
    const md = materialData.material_data || {};
    const nMolecules = md.molecules_hint || 64;
    const density = md.density_g_cm3 || 1.0;

    // Rough box size estimate from density and molecule count
    const boxSize = Math.cbrt(nMolecules / (density * 0.03)) ; // approximate

    const atoms = [];
    for (let i = 0; i < nMolecules; i++) {
      const cx = (Math.random() - 0.5) * boxSize;
      const cy = (Math.random() - 0.5) * boxSize;
      const cz = (Math.random() - 0.5) * boxSize;

      // Place an O with two H's in a rough water-like geometry
      atoms.push({ element: 'O', x: cx, y: cy, z: cz, radius: VDW_RADII.O, color: CPK_COLORS.O });
      atoms.push({ element: 'H', x: cx + 0.76, y: cy + 0.59, z: cz, radius: VDW_RADII.H, color: CPK_COLORS.H });
      atoms.push({ element: 'H', x: cx - 0.76, y: cy + 0.59, z: cz, radius: VDW_RADII.H, color: CPK_COLORS.H });
    }

    return {
      atoms,
      box: { a: boxSize, b: boxSize, c: boxSize },
      metadata: {
        visualization_mode: 'liquid',
        description: materialData.description || 'Liquid phase material',
        atom_count: atoms.length,
        molecule_count: nMolecules
      }
    };
  }

  _buildGas(materialData) {
    const md = materialData.material_data || {};
    const components = md.components || [{ formula: md.molecular_formula || 'N2', fraction: 1.0 }];

    const boxSize = 30; // Angstroms - large sparse box for gas
    const totalMolecules = 32;
    const atoms = [];

    for (const comp of components) {
      const nMol = Math.max(1, Math.round(totalMolecules * (comp.fraction || 1.0)));
      for (let i = 0; i < nMol; i++) {
        const cx = (Math.random() - 0.5) * boxSize;
        const cy = (Math.random() - 0.5) * boxSize;
        const cz = (Math.random() - 0.5) * boxSize;

        const moleculeAtoms = this._gasAtoms(comp.formula, cx, cy, cz);
        atoms.push(...moleculeAtoms);
      }
    }

    return {
      atoms,
      box: { a: boxSize, b: boxSize, c: boxSize },
      metadata: {
        visualization_mode: 'gas',
        description: materialData.description || 'Gas phase material',
        atom_count: atoms.length
      }
    };
  }

  _gasAtoms(formula, cx, cy, cz) {
    const f = (formula || '').toUpperCase();
    if (f === 'N2') {
      return [
        { element: 'N', x: cx - 0.55, y: cy, z: cz, radius: VDW_RADII.N, color: CPK_COLORS.N },
        { element: 'N', x: cx + 0.55, y: cy, z: cz, radius: VDW_RADII.N, color: CPK_COLORS.N }
      ];
    }
    if (f === 'O2') {
      return [
        { element: 'O', x: cx - 0.60, y: cy, z: cz, radius: VDW_RADII.O, color: CPK_COLORS.O },
        { element: 'O', x: cx + 0.60, y: cy, z: cz, radius: VDW_RADII.O, color: CPK_COLORS.O }
      ];
    }
    if (f === 'CO2') {
      return [
        { element: 'C', x: cx, y: cy, z: cz, radius: VDW_RADII.C, color: CPK_COLORS.C },
        { element: 'O', x: cx - 1.16, y: cy, z: cz, radius: VDW_RADII.O, color: CPK_COLORS.O },
        { element: 'O', x: cx + 1.16, y: cy, z: cz, radius: VDW_RADII.O, color: CPK_COLORS.O }
      ];
    }
    if (f === 'H2O') {
      return [
        { element: 'O', x: cx, y: cy, z: cz, radius: VDW_RADII.O, color: CPK_COLORS.O },
        { element: 'H', x: cx + 0.76, y: cy + 0.59, z: cz, radius: VDW_RADII.H, color: CPK_COLORS.H },
        { element: 'H', x: cx - 0.76, y: cy + 0.59, z: cz, radius: VDW_RADII.H, color: CPK_COLORS.H }
      ];
    }
    if (f === 'AR' || f === 'HE' || f === 'NE' || f === 'KR' || f === 'XE') {
      const el = f.charAt(0) + f.slice(1).toLowerCase();
      return [
        { element: el, x: cx, y: cy, z: cz, radius: VDW_RADII[el] || 1.5, color: CPK_COLORS[el] || '#80D1E3' }
      ];
    }
    // Default: single atom placeholder
    return [
      { element: 'C', x: cx, y: cy, z: cz, radius: VDW_RADII.C, color: CPK_COLORS.C }
    ];
  }
}

module.exports = MaterialSceneBuilder;
