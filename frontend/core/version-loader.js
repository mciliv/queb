// Version Loader - Toggle between React and Vanilla JS implementations
class VersionLoader {
  constructor(config = {}) {
    // Configuration options
    this.config = {
      defaultVersion: config.defaultVersion || 'vanilla', // 'vanilla' | 'react'
      allowToggle: config.allowToggle !== false, // true/false
      showToggleButton: config.showToggleButton !== false, // true/false
      persistChoice: config.persistChoice !== false, // true/false
      ...config
    };
    
    this.currentVersion = this.config.persistChoice 
      ? (localStorage.getItem('mol-version') || this.config.defaultVersion)
      : this.config.defaultVersion;
    
    this.isLoading = false;
    this.init();
  }

  init() {
    // Load the default version (React)
    this.loadVersion(this.currentVersion);
  }



  async switchVersion(version) {
    if (!this.config.allowToggle) {
      console.warn('⚠️ Version switching is disabled in configuration');
      return;
    }
    
    if (this.isLoading || version === this.currentVersion) return;
    
    this.isLoading = true;

    try {
      // Clear current app
      this.clearCurrentApp();
      
      // Load new version
      await this.loadVersion(version);
      
      // Update state
      this.currentVersion = version;
      if (this.config.persistChoice) {
        localStorage.setItem('mol-version', version);
      }
      
      console.log(`✅ Switched to ${version} version`);
      
    } catch (error) {
      console.error(`❌ Failed to switch to ${version}:`, error);
      // Revert to previous version
      await this.loadVersion(this.currentVersion);
    } finally {
      this.isLoading = false;
    }
  }

  // Force load a specific version (bypasses toggle restrictions)
  async forceLoadVersion(version) {
    console.log(`🔄 Force loading ${version} version...`);
    
    this.isLoading = true;
    try {
      this.clearCurrentApp();
      await this.loadVersion(version);
      this.currentVersion = version;
      
      if (this.config.persistChoice) {
        localStorage.setItem('mol-version', version);
      }
      
      console.log(`✅ Force loaded ${version} version`);
    } catch (error) {
      console.error(`❌ Failed to force load ${version}:`, error);
      throw error;
    } finally {
      this.isLoading = false;
    }
  }

  // Get current configuration
  getConfig() {
    return { ...this.config };
  }

  // Update configuration at runtime
  updateConfig(newConfig) {
    this.config = { ...this.config, ...newConfig };
    console.log('⚙️ Configuration updated:', this.config);
  }

  clearCurrentApp() {
    const root = document.getElementById('root');
    if (root) {
      root.innerHTML = '';
    }
    
    // Clean up React if it was loaded
    if (window.ReactDOM && window.ReactDOM.unmountComponentAtNode) {
      window.ReactDOM.unmountComponentAtNode(root);
    }
    
    // Clean up global app instance
    if (window.app) {
      if (window.app.cleanup) {
        window.app.cleanup();
      }
      window.app = null;
    }
    
    // Remove any existing event listeners
    document.removeEventListener('imageAnalysisComplete', null);
  }

  async loadVersion(version) {
    console.log(`🔄 Loading ${version} version...`);
    
    if (version === 'react') {
      await this.loadReactVersion();
    } else {
      await this.loadVanillaVersion();
    }
  }

  async loadReactVersion() {
    try {
      // Load React and ReactDOM
      const React = await import('react');
      const { createRoot } = await import('react-dom/client');
      
      // Load React app
      const { default: App } = await import('./App.jsx');
      
      const rootElement = document.getElementById('root');
      const root = createRoot(rootElement);
      const appElement = React.createElement(App);
      
      root.render(appElement);
      
    } catch (error) {
      console.error('Failed to load React version:', error);
      throw error;
    }
  }

  async loadVanillaVersion() {
    try {
      console.log('🔄 Loading Vanilla JS version...');
      
      // Load vanilla JS app
      const { default: MolecularApp } = await import('./app.js');
      console.log('✅ MolecularApp class loaded');
      
      const app = new MolecularApp();
      console.log('✅ MolecularApp instance created');
      
      await app.initialize();
      console.log('✅ Vanilla JS app initialized');
      
      // Make app globally available
      window.app = app;
      
    } catch (error) {
      console.error('Failed to load Vanilla JS version:', error);
      throw error;
    }
  }

  loadScript(src) {
    return new Promise((resolve, reject) => {
      const script = document.createElement('script');
      script.src = src;
      script.onload = resolve;
      script.onerror = reject;
      document.head.appendChild(script);
    });
  }
}

// React-only configuration - no toggle UI
const config = {
  defaultVersion: 'react',
  allowToggle: false,
  showToggleButton: false,
  persistChoice: false
};

// Initialize version loader with configuration
window.versionLoader = new VersionLoader(config);

// Global helper functions
window.switchToReact = () => {
  const loader = window.versionLoader;
  if (loader) loader.switchVersion('react');
};

window.switchToVanilla = () => {
  const loader = window.versionLoader;
  if (loader) loader.switchVersion('vanilla');
};

window.getCurrentVersion = () => {
  return localStorage.getItem('mol-version') || 'vanilla';
};

// Configuration functions
window.setVersionConfig = (config) => {
  if (window.versionLoader) {
    window.versionLoader.updateConfig(config);
  }
  window.VERSION_CONFIG = { ...window.VERSION_CONFIG, ...config };
};

window.getVersionConfig = () => {
  return window.VERSION_CONFIG;
};

window.forceLoadVersion = (version) => {
  const loader = window.versionLoader;
  if (loader) return loader.forceLoadVersion(version);
};

// Test functions for development
window.testVersionSwitch = async () => {
  console.log('🧪 Testing version switch...');
  const currentVersion = window.getCurrentVersion();
  const newVersion = currentVersion === 'react' ? 'vanilla' : 'react';
  
  console.log(`Current: ${currentVersion}, Switching to: ${newVersion}`);
  await window.versionLoader.switchVersion(newVersion);
  
  console.log('✅ Version switch test completed');
};

window.logVersionInfo = () => {
  console.log('📊 Version Information:');
  console.log('- Current version:', window.getCurrentVersion());
  console.log('- Version loader:', window.versionLoader ? 'Loaded' : 'Not loaded');
  console.log('- Configuration:', window.getVersionConfig());
  console.log('- React available:', typeof window.React !== 'undefined');
  console.log('- Vanilla app:', window.app ? 'Loaded' : 'Not loaded');
}; 