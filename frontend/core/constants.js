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

export const PRESET_VISUAL_TESTS = [
  { label: TEST_MOLECULES.water.name, smiles: TEST_MOLECULES.water.smiles },
  { label: TEST_MOLECULES.ethanol.name, smiles: TEST_MOLECULES.ethanol.smiles },
  { label: TEST_MOLECULES.caffeine.name, smiles: TEST_MOLECULES.caffeine.smiles }
];
