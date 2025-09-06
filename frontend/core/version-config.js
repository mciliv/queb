// Client-side version configuration
// No server-side dependencies for browser compatibility

export const VERSION_CONFIG = {
  defaultVersion: 'react',
  allowToggle: false,
  showToggleButton: false,
  persistChoice: false,
  performanceMode: false
};

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