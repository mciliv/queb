// Client-side constants for browser compatibility
// Loads configuration from external config file

// Use hardcoded fallback data for better compatibility
const configData = {
  TEST_MOLECULES: {
    water: { name: 'Water', smiles: 'O' },
    ethanol: { name: 'Ethanol', smiles: 'CCO' },
    caffeine: { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' }
  },
  SMILES_NAME_MAP: {
    'O': 'Water',
    'CCO': 'Ethanol',
    'CN1C=NC2=C1C(=O)N(C(=O)N2C)C': 'Caffeine'
  },
  PRESET_VISUAL_TESTS: [
    {
      label: 'Coffee',
      smilesList: ['CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'O', 'CCO']
    },
    {
      label: 'Red wine',
      smilesList: ['O', 'CCO', 'CC(=O)O', 'C1=CC=CC=C1']
    },
    {
      label: 'Apple',
      smilesList: ['O', 'CC(=O)O']
    }
  ]
};

export const TEST_MOLECULES = configData.TEST_MOLECULES;
export const SMILES_NAME_MAP = configData.SMILES_NAME_MAP;
export const PRESET_VISUAL_TESTS = configData.PRESET_VISUAL_TESTS;
