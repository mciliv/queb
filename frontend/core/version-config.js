// Version Configuration - now consolidated in project.js
// This file is deprecated, use project.js configuration instead

import project from '../../config/project.js';

export const VERSION_CONFIG = project.helpers.getVersionConfig();
export const getActiveConfig = () => VERSION_CONFIG;
export const getFinalConfig = () => VERSION_CONFIG;
export const CONFIG_PRESETS = {};
export const applyPreset = () => false;

// Export for global access
window.VERSION_CONFIG = VERSION_CONFIG;
window.getActiveConfig = getActiveConfig;
window.getFinalConfig = getFinalConfig;
window.CONFIG_PRESETS = CONFIG_PRESETS;
window.applyPreset = applyPreset; 