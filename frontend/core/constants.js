// Client-side constants for browser compatibility
// Loads configuration from external config file

import configData from '../config/molecular-config.json' with { type: 'json' };

export const TEST_MOLECULES = configData.TEST_MOLECULES;
export const SMILES_NAME_MAP = configData.SMILES_NAME_MAP;
export const PRESET_VISUAL_TESTS = configData.PRESET_VISUAL_TESTS;
