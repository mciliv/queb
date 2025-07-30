// Simple App Logic - Following ui.mdc Guidelines
import { simplePaymentManager } from '../components/simple-payment.js';
import { cameraManager } from '../components/camera.js';
import { cameraHandler } from '../components/camera-handler.js';
import { uiManager } from '../components/ui-utils.js';
import { logger } from '../components/logger.js';

class MolecularApp {
  constructor() {
    this.snapshots = null;
    this.objectInput = null;
    this.viewers = [];
    this.isProcessing = false;
    this.hasPaymentSetup = false;
    this.currentAnalysisType = null;
    this.lastAnalysis = null;
  }

  async initialize() {
    this.snapshots = document.querySelector(".snapshots-container");
    this.objectInput = document.getElementById("object-input");

    uiManager.initialize();
    uiManager.setupDebuggingFunctions();
    uiManager.showMainApp();

    this.setupEventListeners();

    // Initialize sidebar payment system
    await simplePaymentManager.checkPaymentRequired();
    
    // Auto-enable dev mode for localhost
    if (location.hostname === 'localhost' || location.hostname === '127.0.0.1') {
      logger.info('Auto-enabling developer mode for localhost');
      this.hasPaymentSetup = true;
    }

    logger.info('Molecular analysis app initialized');
  }

  setupEventListeners() {
    this.setupTextAnalysis();
    cameraHandler.setupEventListeners();
    
    // Keyboard shortcut (Cmd+K / Ctrl+K)
    document.addEventListener('keydown', (event) => {
      if ((event.metaKey || event.ctrlKey) && event.key === 'k') {
        event.preventDefault();
        const textInput = document.getElementById('object-input');
        if (textInput) {
          textInput.focus();
          textInput.select();
        }
      }
    });
    
    // Platform-specific placeholder
    const textInput = document.getElementById('object-input');
    if (textInput) {
      const isMac = navigator.platform.toUpperCase().indexOf('MAC') >= 0;
      const shortcutKey = isMac ? '⌘K' : 'Ctrl+K';
      textInput.placeholder = `Describe an object... (${shortcutKey} to focus)`;
    }

    document.addEventListener('imageAnalysisComplete', (e) => {
      const { output, icon, objectName, useQuotes, croppedImageData } = e.detail;
      this.processAnalysisResult(output, icon, objectName, useQuotes, croppedImageData);
    });

    // Camera mode selection handler
    const cameraMode = document.getElementById('camera-mode');
    if (cameraMode) {
      cameraMode.addEventListener('change', async (e) => {
        if (e.target.checked) {
          logger.cameraEvent('camera_mode_activated');
          try {
            await cameraManager.initialize();
            const permissionGranted = await cameraManager.requestPermission();
            if (!permissionGranted) {
              e.target.checked = false;
            }
          } catch (error) {
            logger.error('Camera initialization failed', error);
            e.target.checked = false;
          }
        }
      });
    }

    // Card management button handler
    const cardBtn = document.getElementById('card-icon-btn');
    if (cardBtn) {
      cardBtn.addEventListener('click', () => {
        logger.userAction('payment_sidebar_toggle');
        const paymentSection = document.getElementById('payment-section');
        if (paymentSection) {
          paymentSection.classList.toggle('collapsed');
        }
      });
    }
  }

  setupTextAnalysis() {
    const textInput = document.getElementById("object-input");
    if (textInput) {
      // Prevent multiple rapid triggers with debouncing
      let lastTriggerTime = 0;
      const DEBOUNCE_MS = 1000; // 1 second cooldown
      
      textInput.addEventListener("keydown", (event) => {
        if (event.key === "Enter") {
          event.preventDefault();
          
          const now = Date.now();
          if (now - lastTriggerTime < DEBOUNCE_MS) {
            
            return;
          }
          
          lastTriggerTime = now;
          this.handleTextAnalysis();
        }
      });
    }
  }
  async handleTextAnalysis() {
    if (this.isProcessing) return;

    const inputValue = this.objectInput.value.trim();
    if (!inputValue) return;

    // Check payment setup for sidebar
    const paymentSetup = await this.checkPaymentSetupForSidebar();
    if (!paymentSetup) {
      logger.warn('Payment not set up, showing message');
      this.showError('Payment setup required. Complete setup in the sidebar on the right.');
      return;
    }

    this.isProcessing = true;
    this.currentAnalysisType = 'text';
    
    try {
      this.showProcessing();
      
      logger.analysisEvent('text_analysis_started', { input: inputValue });
      
      // Create AbortController for timeout
      const controller = new AbortController();
      const timeoutId = setTimeout(() => controller.abort(), 45000); // 45 second timeout
      
      const response = await fetch("/analyze-text", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ object: inputValue }),
        signal: controller.signal
      });

      clearTimeout(timeoutId);

      if (!response.ok) {
        const errorText = await response.text();
        throw new Error(`Failed to analyze text: ${response.status} ${errorText}`);
      }

      const result = await response.json();
      logger.analysisEvent('text_analysis_completed', { input: inputValue, result });

      this.lastAnalysis = {
        type: 'text',
        input: inputValue,
        result: result
      };

      this.processAnalysisResult(result.output, null, inputValue, false, null);
      
    } catch (error) {
      logger.error('Text analysis failed', error);
      
      // Enhanced error handling with user-friendly messages
      let errorMessage = 'Analysis failed: ';
      
      if (error.name === 'AbortError') {
        errorMessage = 'Analysis timed out. Please check your internet connection and try again.';
      } else if (error.message.includes('Failed to fetch') || error.message.includes('fetch')) {
        errorMessage = 'Network connection failed. Please check your internet connection and try again.';
      } else if (error.message.includes('Network connection failed')) {
        errorMessage = 'Unable to connect to analysis service. Please check your internet connection.';
      } else if (error.message.includes('Rate limit exceeded')) {
        errorMessage = 'Too many requests. Please wait a moment and try again.';
      } else if (error.message.includes('503') || error.message.includes('temporarily unavailable')) {
        errorMessage = 'Analysis service temporarily unavailable. Please try again in a few moments.';
      } else {
        errorMessage += error.message;
      }
      
      this.showError(errorMessage);
    } finally {
      this.hideProcessing();
      this.isProcessing = false;
    }
  }

  // Check payment setup for sidebar-based system
  async checkPaymentSetupForSidebar() {
    // Check if dev mode is enabled (localhost auto-enable)
    if (this.hasPaymentSetup === true) {
      logger.info('Developer mode active - bypassing payment check');
      return true;
    }
    
    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    if (!deviceToken || !cardInfo) {
      // Ensure payment section is visible
      const paymentSection = document.getElementById('payment-section');
      if (paymentSection) {
        paymentSection.classList.remove('hidden');
      }
      return false;
    }
    
    try {
      const response = await fetch('/validate-payment', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ device_token: deviceToken })
      });
      
      if (!response.ok) {
        localStorage.removeItem('molDeviceToken');
        localStorage.removeItem('molCardInfo');
        // Ensure payment section is visible
        const paymentSection = document.getElementById('payment-section');
        if (paymentSection) {
          paymentSection.classList.remove('hidden');
        }
        return false;
      }
      
      return true;
      
    } catch (error) {
      logger.error('Payment validation error', error);
      return true; // Fallback to allow analysis
    }
  }

  showProcessing() {
    // Create loading indicator
    const loadingDiv = document.createElement('div');
    loadingDiv.className = 'loading-indicator';
    loadingDiv.innerHTML = `
      <div class="spinner"></div>
      <div>Analyzing...</div>
    `;
    loadingDiv.style.cssText = `
      position: fixed;
      top: 50%;
      left: 50%;
      transform: translate(-50%, -50%);
      background: rgba(0, 0, 0, 0.9);
      color: white;
      padding: 20px;
      border-radius: 10px;
      text-align: center;
      z-index: 1000;
    `;
    
    document.body.appendChild(loadingDiv);
    this.loadingIndicator = loadingDiv;
  }

  hideProcessing() {
    if (this.loadingIndicator) {
      this.loadingIndicator.remove();
      this.loadingIndicator = null;
    }
  }

  showError(message) {
    logger.error('Application error', { message });
    // Create simple error display
    const errorDiv = document.createElement('div');
    errorDiv.className = 'error-message';
    errorDiv.textContent = message;
    errorDiv.style.cssText = `
      background: rgba(255, 68, 68, 0.1);
      border: 1px solid #ff4444;
      border-radius: 4px;
      padding: 12px;
      margin: 12px 16px;
      color: #ff4444;
      font-size: 12px;
      text-align: center;
    `;
    
    const resultsSection = document.querySelector('.results-section');
    if (resultsSection) {
      // Remove existing error messages
      const existingErrors = resultsSection.querySelectorAll('.error-message');
      existingErrors.forEach(err => err.remove());
      
      resultsSection.insertBefore(errorDiv, resultsSection.firstChild);
      
      // Auto-remove after 5 seconds
      setTimeout(() => errorDiv.remove(), 5000);
    }
  }

  // Simple analysis result processing
  async processAnalysisResult(output, icon, objectName, useQuotes = false, croppedImageData = null) {
    const chemicals = output.chemicals || [];
    const displayName = useQuotes ? `"${objectName}"` : objectName;
    
    // Extract valid SMILES
    const smiles = chemicals
      .map(c => c.smiles)
      .filter(s => s && s !== "N/A" && s.trim() !== "");
    
    // Generate and display molecules
    if (smiles.length > 0) {
      await this.generateAndDisplayMolecules(smiles, displayName, chemicals);
    } else {
      // Better error message based on what was found
      let message = "No displayable molecular structures found";
      if (chemicals.length > 0) {
        message = `Found ${chemicals.length} chemical(s) but no valid SMILES structures for 3D visualization`;
      }
      await this.createObjectColumn(displayName, [], [], message);
    }
  }

  // Generate SDF files and create molecular display
  async generateAndDisplayMolecules(smiles, objectName, chemicals) {
    try {
      const response = await fetch("/generate-sdfs", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ smiles, overwrite: true }),
      });

      if (!response.ok) throw new Error(`SDF generation failed: ${response.status}`);
      
      const result = await response.json();
      
      // Check if any molecules were successfully generated
      if (result.sdfPaths && result.sdfPaths.length > 0) {
        await this.createObjectColumn(objectName, result.sdfPaths, smiles, null, chemicals);
      } else {
        // Show what went wrong
        let message = "No 3D structures could be generated";
        if (result.skipped && result.skipped.length > 0) {
          message += `\n\nSkipped: ${result.skipped.join(', ')}`;
        }
        if (result.errors && result.errors.length > 0) {
          message += `\n\nErrors: ${result.errors.join(', ')}`;
        }
        await this.createObjectColumn(objectName, [], [], message);
      }
      
    } catch (error) {
      this.showError(`Error generating 3D models: ${error.message}`);
    }
  }

  // Create molecular display column
  async createObjectColumn(objectName, sdfFiles, smiles = [], description = null, chemicals = null) {
    const gldiv = document.getElementById("gldiv");
    if (!gldiv) return;

    const objectColumn = document.createElement("div");
    objectColumn.className = "object-column";

    // Title with close button
    const titleContainer = document.createElement("div");
    titleContainer.className = "object-title";
    titleContainer.innerHTML = `
      <span>${objectName}</span>
      <button class="close-button" onclick="this.parentElement.parentElement.remove()">✕</button>
    `;
    objectColumn.appendChild(titleContainer);

    if (description) {
      // Show description
      const descDiv = document.createElement("div");
      descDiv.className = "description-content";
      descDiv.textContent = description;
      objectColumn.appendChild(descDiv);
    } else {
      // Render 3D molecules
      for (let i = 0; i < sdfFiles.length; i++) {
        const moleculeContainer = document.createElement("div");
        moleculeContainer.className = "molecule-container";
        
        const moleculeName = document.createElement("div");
        moleculeName.className = "molecule-name";
        const chemical = chemicals?.[i] || { smiles: smiles[i] };
        moleculeName.textContent = uiManager.getMoleculeName(chemical);
        
        const viewerContainer = document.createElement("div");
        viewerContainer.className = "molecule-viewer";
        viewerContainer.id = `viewer-${Date.now()}-${i}`;
        
        moleculeContainer.appendChild(moleculeName);
        moleculeContainer.appendChild(viewerContainer);
        objectColumn.appendChild(moleculeContainer);

        // Add small delay to ensure DOM is ready
        setTimeout(async () => {
          const viewer = await this.render(sdfFiles[i], viewerContainer);
          if (viewer) this.viewers.push(viewer);
        }, i * 100); // Stagger rendering to avoid conflicts
      }
    }

    gldiv.appendChild(objectColumn);
  }

  async render(sdfFile, container) {
    try {
      const response = await fetch(sdfFile);
      if (!response.ok) throw new Error(`HTTP error ${response.status}`);

      const sdfData = await response.text();
      
      // Clear container first
      container.innerHTML = '';
      
      // Create viewer with proper configuration
      const viewer = $3Dmol.createViewer(container, { 
        defaultcolors: $3Dmol.rasmolElementColors,
        backgroundColor: 'black'
      });
      
      viewer.addModel(sdfData, "sdf");
      viewer.setStyle({}, { sphere: { scale: 0.8 } });
      viewer.zoomTo();
      viewer.render();
      
      // Ensure viewer resizes properly
      setTimeout(() => {
        viewer.resize();
        viewer.render();
      }, 50);
      
      return viewer;
    } catch (error) {
      container.innerHTML = `<div style="color: #ff4444; padding: 10px; text-align: center;">Error: ${error.message}</div>`;
      return null;
    }
  }

  // Alias for test compatibility
  handleError(error) {
    this.showError(error.message);
  }

  // Generate SDFs (alias for generateAndDisplayMolecules for test compatibility)
  async generateSDFs(smiles, objectName, icon, chemicals, croppedImageData) {
    try {
      const response = await fetch("/generate-sdfs", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ smiles, overwrite: true }),
      });

      if (!response.ok) throw new Error(`SDF generation failed: ${response.status}`);
      
      const result = await response.json();
      await this.createObjectColumn(objectName, result.sdfPaths || [], smiles, null, chemicals);
      
    } catch (error) {
      this.showError(`Error generating 3D models: ${error.message}`);
    }
  }

  // Create closable error message (for test compatibility)
  createClosableErrorMessage(message) {
    return uiManager.createErrorMessage(message, this.snapshots);
  }

  // Cleanup method for test compatibility
  cleanup() {
    // Clean up viewers
    this.viewers.forEach(viewer => {
      if (viewer && typeof viewer.clear === 'function') {
        viewer.clear();
      }
    });
    this.viewers = [];

    // Clean up managers if they have cleanup methods
    if (cameraManager && typeof cameraManager.cleanup === 'function') {
      cameraManager.cleanup();
    }
    if (uiManager && typeof uiManager.cleanup === 'function') {
      uiManager.cleanup();
    }
  }

  // Processing indicator methods for test compatibility
  showProcessing() {
    const indicator = document.getElementById('processing-indicator');
    if (indicator) indicator.style.display = 'block';
  }

  hideProcessing() {
    const indicator = document.getElementById('processing-indicator');
    if (indicator) indicator.style.display = 'none';
  }

  updateScrollHandles() {
    // Simple implementation for horizontal scroll management
    const gldiv = document.getElementById("gldiv");
    if (gldiv) {
      const hasOverflow = gldiv.scrollWidth > gldiv.clientWidth;
      gldiv.style.overflowX = hasOverflow ? 'auto' : 'hidden';
    }
  }

  // Debug layout toggle (development only)
  enableLayoutDebug() {
    document.body.classList.toggle('debug-layout');
     ? 'ON' : 'OFF');
  }

  // UI Context Automation integration (manual activation only)
  initUIContextAutomation() {
    if (window.uiContextAutomation) {
       to activate manually.');
    }
  }

}

// Export for testing
export { MolecularApp };

// Initialize app when DOM is ready
document.addEventListener("DOMContentLoaded", async () => {
  const app = new MolecularApp();
  window.app = app; // Make app globally available for debugging
  await app.initialize();

  window.molecularApp = app;
  
  // Initialize UI Context Automation after app is created (manual activation only)
  if (window.uiContextAutomation) {
    window.app.enableUIContextAutomation = () => window.uiContextAutomation.enable();
    window.app.captureUIState = (label) => window.uiContextAutomation.captureUIState(label);
    window.app.exportContextForPrompt = () => window.uiContextAutomation.exportContextForPrompt();
    window.app.checkRuleCompliance = () => window.uiContextAutomation.checkRuleCompliance();
    
    // UI Context Automation available but not auto-enabled (prevents Chrome crashes)
     to activate manually.');
  }
});

