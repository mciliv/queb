// Frontend constants and test data
export const TEST_MOLECULES = {
  water: {
    name: 'Water',
    smiles: 'O'
  },
  ethanol: {
    name: 'Ethanol', 
    smiles: 'CCO'
  },
  caffeine: {
    name: 'Caffeine',
    smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'
  }
};

// Visual tests should be objects with many compounds (use multiple SMILES per object)
// Use SMILES we can reliably generate offline from existing fixtures
export const PRESET_VISUAL_TESTS = [
  {
    label: 'Coffee',
    smilesList: [
      TEST_MOLECULES.caffeine.smiles, // caffeine
      TEST_MOLECULES.water.smiles,    // water
      TEST_MOLECULES.ethanol.smiles   // trace ethanol (for test coverage)
    ]
  },
  {
    label: 'Red wine',
    smilesList: [
      TEST_MOLECULES.water.smiles,
      TEST_MOLECULES.ethanol.smiles,
      'CC(=O)O',            // acetic acid
      'C1=CC=CC=C1'         // benzene (phenolic proxy for test)
    ]
  },
  {
    label: 'Apple',
    smilesList: [
      TEST_MOLECULES.water.smiles,
      'CC(=O)O'             // acetic acid (vinegar note) as placeholder
    ]
  }
];
