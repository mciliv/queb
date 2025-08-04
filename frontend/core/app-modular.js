// app-modular.js - Modular application with component-based architecture
import { AppShell } from './components/AppShell.js';
import { MolecularViewer } from './components/MolecularViewer.js';
import { ErrorHandler } from './components/ErrorHandler.js';
import { simplePaymentManager } from '../components/simple-payment.js';
import { cameraManager } from '../components/camera.js';
import { cameraHandler } from '../components/camera-handler.js';
import { uiManager } from '../components/ui-utils.js';

// Payment toggle configuration
const PAYMENT_CONFIG = {
  enabled: false, // Set to true to enable payment functionality
  devMode: location.hostname === 'localhost' || location.hostname === '127.0.0.1'
};

class ModularMolecularApp {
  constructor() {
    this.appShell = new AppShell();
    this.molecularViewer = new MolecularViewer();
    this.errorHandler = new ErrorHandler();
    this.paymentEnabled = PAYMENT_CONFIG.enabled;
    this.hasPaymentSetup = false;
  }

  async initialize() {
    try {
      // Initialize UI manager first
      uiManager.initialize();
      uiManager.setupDebuggingFunctions();
      uiManager.showMainApp();

      // Initialize error handler
      this.errorHandler.initialize();
      this.errorHandler.setupGlobalErrorHandler();

      // Initialize app shell (handles keyboard shortcuts)
      await this.appShell.initialize();

      // Initialize molecular viewer
      await this.molecularViewer.initialize();

      // Initialize payment system if enabled
      if (this.paymentEnabled) {
        simplePaymentManager.checkPaymentRequired();
        this.hasPaymentSetup = true;
      } else {
        // Hide payment section when disabled
        this.appShell.hidePaymentSection();
        this.hasPaymentSetup = true; // Skip payment validation
        console.log('💳 Payment functionality disabled');
      }
      
      // Auto-enable dev mode for localhost
      if (PAYMENT_CONFIG.devMode) {
        console.log('🔧 Auto-enabling developer mode for localhost');
        this.hasPaymentSetup = true;
      }    

      // Initialize camera system
      await cameraManager.initialize();
      cameraHandler.setupEventListeners();

      if (cameraManager.isSafari && !cameraManager.hasStoredCameraPermission()) {
        setTimeout(() => {
          cameraManager.requestPermission();
        }, 1000);
      }

      console.log('✅ Modular molecular analysis app initialized');
    } catch (error) {
      console.error('❌ Failed to initialize app:', error);
      this.errorHandler.handleError('Failed to initialize application', 'critical');
    }
  }

  // Public methods for external access
  getAppShell() {
    return this.appShell;
  }

  getMolecularViewer() {
    return this.molecularViewer;
  }

  getErrorHandler() {
    return this.errorHandler;
  }

  // Payment management
  togglePayments(enabled) {
    this.paymentEnabled = enabled;
    this.appShell.setPaymentStatus(enabled, this.hasPaymentSetup);
    
    if (enabled) {
      this.appShell.showPaymentSection();
      simplePaymentManager.checkPaymentRequired();
    } else {
      this.appShell.hidePaymentSection();
    }
  }

  // Clear all results
  clearResults() {
    this.appShell.clearResults();
    this.molecularViewer.clearResults();
  }
}

// Initialize the app when DOM is ready
document.addEventListener('DOMContentLoaded', async () => {
  window.molecularApp = new ModularMolecularApp();
  await window.molecularApp.initialize();
});

// Export for module usage
export { ModularMolecularApp }; 