// Client-side constants for browser compatibility
// Loads configuration from external config file

import fs from 'fs';
import path from 'path';

// Load molecular configuration from JSON file
const configPath = path.join(process.cwd(), 'frontend', 'config', 'molecular-config.json');
const configData = JSON.parse(fs.readFileSync(configPath, 'utf8'));

export const TEST_MOLECULES = configData.TEST_MOLECULES;
export const SMILES_NAME_MAP = configData.SMILES_NAME_MAP;
export const PRESET_VISUAL_TESTS = configData.PRESET_VISUAL_TESTS;
