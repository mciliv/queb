// Version Loader - Toggle between React and Vanilla JS implementations
class VersionLoader {
  constructor(config = {}) {
    console.log('🔧 VersionLoader constructor called with config:', config);
    
    // Configuration options
    this.config = {
      defaultVersion: config.defaultVersion || 'vanilla', // 'vanilla' | 'react'
      allowToggle: config.allowToggle !== false, // true/false
      showToggleButton: config.showToggleButton !== false, // true/false
      persistChoice: config.persistChoice !== false, // true/false
      ...config
    };
    
    console.log('🔧 Final config:', this.config);
    
    this.currentVersion = this.config.persistChoice 
      ? (localStorage.getItem('mol-version') || this.config.defaultVersion)
      : this.config.defaultVersion;
    
    console.log('🔧 Current version will be:', this.currentVersion);
    this.isLoading = false;
    this.init();
  }

  init() {
    if (this.config.showToggleButton) {
      this.setupToggle();
    }
    // No toggle button setup needed - React is the default
    this.loadVersion(this.currentVersion);
  }

  setupToggle() {
    const toggle = document.getElementById('version-toggle');
    if (!toggle) return;

    // Only add event listener if toggling is allowed
    if (this.config.allowToggle) {
      toggle.addEventListener('click', () => {
        if (this.isLoading) return;
        
        const newVersion = this.currentVersion === 'react' ? 'vanilla' : 'react';
        this.switchVersion(newVersion);
      });
    } else {
      toggle.style.cursor = 'not-allowed';
      toggle.title = 'Version switching disabled';
    }

    // Update toggle appearance
    this.updateToggleAppearance();
  }

  updateToggleAppearance() {
    const toggle = document.getElementById('version-toggle');
    if (!toggle) return; // No toggle button to update

    toggle.className = `version-toggle ${this.currentVersion}`;
    toggle.textContent = this.currentVersion === 'react' ? 'React' : 'Vanilla JS';
    
    if (this.config.allowToggle) {
      toggle.title = `Switch to ${this.currentVersion === 'react' ? 'Vanilla JS' : 'React'}`;
    } else {
      toggle.title = 'Version switching disabled';
    }
  }

  async switchVersion(version) {
    if (!this.config.allowToggle) {
      console.warn('⚠️ Version switching is disabled in configuration');
      return;
    }
    
    if (this.isLoading || version === this.currentVersion) return;
    
    this.isLoading = true;
    const toggle = document.getElementById('version-toggle');
    if (toggle) {
      toggle.textContent = 'Switching...';
      toggle.style.cursor = 'not-allowed';
    }

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
      this.updateToggleAppearance();
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
      this.updateToggleAppearance();
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
      console.log('🔄 Loading React version...');
      
      // Load React and ReactDOM
      const React = await import('react');
      const { createRoot } = await import('react-dom/client');
      console.log('✅ React and ReactDOM loaded');
      
      // Load React app
      const { default: App } = await import('./App.jsx');
      console.log('✅ React App component loaded');
      
      const rootElement = document.getElementById('root');
      console.log('✅ Root element found:', rootElement);
      
      const root = createRoot(rootElement);
      console.log('✅ React root created');
      
      const appElement = React.createElement(App);
      console.log('✅ App element created');
      
      root.render(appElement);
      console.log('✅ React app rendered');
      
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

// Simple React default configuration
const config = {
  defaultVersion: 'react',
  allowToggle: false,
  showToggleButton: false,
  persistChoice: false
};

console.log('🔧 Using React default config:', config);

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