// Version Loader - Toggle between React and Vanilla JS implementations
class VersionLoader {
  constructor() {
    this.currentVersion = localStorage.getItem('mol-version') || 'vanilla';
    this.isLoading = false;
    this.init();
  }

  init() {
    this.setupToggle();
    this.loadVersion(this.currentVersion);
  }

  setupToggle() {
    const toggle = document.getElementById('version-toggle');
    if (!toggle) return;

    toggle.addEventListener('click', () => {
      if (this.isLoading) return;
      
      const newVersion = this.currentVersion === 'react' ? 'vanilla' : 'react';
      this.switchVersion(newVersion);
    });

    // Update toggle appearance
    this.updateToggleAppearance();
  }

  updateToggleAppearance() {
    const toggle = document.getElementById('version-toggle');
    if (!toggle) return;

    toggle.className = `version-toggle ${this.currentVersion}`;
    toggle.textContent = this.currentVersion === 'react' ? 'React' : 'Vanilla JS';
    toggle.title = `Switch to ${this.currentVersion === 'react' ? 'Vanilla JS' : 'React'}`;
  }

  async switchVersion(version) {
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
      localStorage.setItem('mol-version', version);
      
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
      
      // Load React app directly (Vite handles React dependencies)
      const { default: App } = await import('./App.jsx');
      console.log('✅ React App component loaded');
      
      const { createRoot } = await import('react-dom/client');
      console.log('✅ ReactDOM createRoot loaded');
      
      const root = createRoot(document.getElementById('root'));
      root.render(App());
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

// Initialize version loader
window.versionLoader = new VersionLoader();

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
  console.log('- React available:', typeof window.React !== 'undefined');
  console.log('- Vanilla app:', window.app ? 'Loaded' : 'Not loaded');
}; 