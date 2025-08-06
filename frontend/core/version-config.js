// Version Configuration - Control React vs Vanilla JS behavior
// Modify these settings to control which version loads and how

export const VERSION_CONFIG = {
  // Default version to load
  defaultVersion: 'vanilla', // 'vanilla' | 'react'
  
  // Toggle functionality
  allowToggle: true,         // Enable/disable version switching
  showToggleButton: true,    // Show/hide the toggle button
  
  // Persistence
  persistChoice: true,       // Remember user's choice in localStorage
  
  // Performance mode (for production)
  performanceMode: false,    // When true, forces vanilla JS, hides toggle
  
  // Environment-specific overrides
  environment: {
    development: {
      allowToggle: true,
      showToggleButton: true,
      defaultVersion: 'vanilla'
    },
    production: {
      allowToggle: false,
      showToggleButton: false,
      defaultVersion: 'vanilla',
      performanceMode: true
    },
    testing: {
      allowToggle: true,
      showToggleButton: true,
      defaultVersion: 'react'
    }
  }
};

// Auto-detect environment and apply overrides
const getEnvironmentConfig = () => {
  const isDev = window.location.hostname === 'localhost' || 
                window.location.hostname === '127.0.0.1' ||
                window.location.hostname.includes('dev');
  
  const isTest = window.location.search.includes('test=true') ||
                 window.location.hostname.includes('test');
  
  if (isTest) {
    return VERSION_CONFIG.environment.testing;
  } else if (isDev) {
    return VERSION_CONFIG.environment.development;
  } else {
    return VERSION_CONFIG.environment.production;
  }
};

// Merge base config with environment-specific config
export const getActiveConfig = () => {
  const envConfig = getEnvironmentConfig();
  return {
    ...VERSION_CONFIG,
    ...envConfig,
    // Remove environment object from final config
    environment: undefined
  };
};

// URL parameter overrides
export const getUrlOverrides = () => {
  const params = new URLSearchParams(window.location.search);
  const overrides = {};
  
  if (params.has('version')) {
    overrides.defaultVersion = params.get('version');
  }
  
  if (params.has('toggle')) {
    overrides.allowToggle = params.get('toggle') === 'true';
  }
  
  if (params.has('showToggle')) {
    overrides.showToggleButton = params.get('showToggle') === 'true';
  }
  
  if (params.has('persist')) {
    overrides.persistChoice = params.get('persist') === 'true';
  }
  
  return overrides;
};

// Get final configuration with all overrides applied
export const getFinalConfig = () => {
  const baseConfig = getActiveConfig();
  const urlOverrides = getUrlOverrides();
  
  return {
    ...baseConfig,
    ...urlOverrides
  };
};

// Configuration presets for common scenarios
export const CONFIG_PRESETS = {
  // Development with full toggle functionality
  development: {
    defaultVersion: 'vanilla',
    allowToggle: true,
    showToggleButton: true,
    persistChoice: true,
    performanceMode: false
  },
  
  // Production optimized for performance
  production: {
    defaultVersion: 'vanilla',
    allowToggle: false,
    showToggleButton: false,
    persistChoice: false,
    performanceMode: true
  },
  
  // Testing with React as default
  testing: {
    defaultVersion: 'react',
    allowToggle: true,
    showToggleButton: true,
    persistChoice: false,
    performanceMode: false
  },
  
  // React-only mode
  reactOnly: {
    defaultVersion: 'react',
    allowToggle: false,
    showToggleButton: false,
    persistChoice: false,
    performanceMode: false
  },
  
  // Vanilla-only mode
  vanillaOnly: {
    defaultVersion: 'vanilla',
    allowToggle: false,
    showToggleButton: false,
    persistChoice: false,
    performanceMode: true
  }
};

// Apply a preset configuration
export const applyPreset = (presetName) => {
  const preset = CONFIG_PRESETS[presetName];
  if (!preset) {
    console.error(`❌ Unknown preset: ${presetName}`);
    return false;
  }
  
  Object.assign(VERSION_CONFIG, preset);
  console.log(`✅ Applied preset: ${presetName}`);
  return true;
};

// Export for global access
window.VERSION_CONFIG = VERSION_CONFIG;
window.getActiveConfig = getActiveConfig;
window.getFinalConfig = getFinalConfig;
window.CONFIG_PRESETS = CONFIG_PRESETS;
window.applyPreset = applyPreset; 