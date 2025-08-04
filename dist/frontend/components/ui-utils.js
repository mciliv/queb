// ui-utils.js - UI utilities and helper functions

class UIManager {
  constructor() {
    this.cameraMode = null;
    this.photoMode = null;
    this.cameraContainer = null;
    this.photoOptions = null;
  }

  // Initialize UI elements
  initialize() {
    this.cameraMode = document.getElementById("camera-mode");
    this.photoMode = document.getElementById("photo-mode");
    this.cameraContainer = document.querySelector(".camera-container");
    this.photoOptions = document.getElementById("photo-options");

    this.setupEventListeners();
    this.updateInputMode();
  }

  // Setup UI event listeners
  setupEventListeners() {
    this.cameraMode.addEventListener("change", () => this.updateInputMode());
    this.photoMode.addEventListener("change", () => this.updateInputMode());
  }

  // Update input mode display
  updateInputMode() {
    // Show/hide camera based on checkbox state
    if (this.cameraMode.checked) {
      this.cameraContainer.classList.add('active');
    } else {
      this.cameraContainer.classList.remove('active');
    }

    // Show/hide photo options based on checkbox state
    if (this.photoMode.checked) {
      this.photoOptions.classList.add('active');
    } else {
      this.photoOptions.classList.remove('active');
    }
  }

  // Mode switching helpers
  switchToCameraMode() {
    this.cameraMode.checked = true;
    this.photoMode.checked = false;
    this.updateInputMode();
  }

  switchToPhotoMode() {
    this.photoMode.checked = true;
    this.cameraMode.checked = false;
    this.updateInputMode();
  }

  clearModeSelection() {
    this.cameraMode.checked = false;
    this.photoMode.checked = false;
    this.updateInputMode();
  }

  // Create error message
  createErrorMessage(message, container) {
    const errorDiv = document.createElement("div");
    errorDiv.className = "error-message-container";

    errorDiv.innerHTML = `
      <div class="error-message-title">⚠️ Error</div>
      <div>${message}</div>
      <button class="error-close-btn" onclick="this.parentElement.remove(); window.updateScrollHandles && window.updateScrollHandles();">×</button>
    `;

    if (container) {
      container.appendChild(errorDiv);
    }

    return errorDiv;
  }

  // Create closable error message
  createClosableErrorMessage(message) {
    const gldiv = document.getElementById("gldiv");
    if (!gldiv) return null;

    const errorColumn = document.createElement("div");
    errorColumn.className = "object-column error-column";

    const header = document.createElement("div");
    header.className = "object-header";
    header.innerHTML = `
      <div class="object-title-container">
        <div class="object-icon">⚠️</div>
        <div class="object-name">Error</div>
      </div>
      <button class="object-close" onclick="this.parentElement.parentElement.remove(); window.updateScrollHandles && window.updateScrollHandles();">
        <img src="close.svg" alt="Close" width="16" height="16" />
      </button>
    `;

    const errorContent = document.createElement("div");
    errorContent.className = "error-content";
    errorContent.textContent = message;

    errorColumn.appendChild(header);
    errorColumn.appendChild(errorContent);
    gldiv.appendChild(errorColumn);

    if (window.updateScrollHandles) {
      window.updateScrollHandles();
    }

    return errorColumn;
  }

  // Create a simple column
  createColumn(message, className = "") {
    const gldiv = document.getElementById("gldiv");
    if (!gldiv) return null;

    const column = document.createElement("div");
    column.className = `object-column ${className}`;

    const header = document.createElement("div");
    header.className = "object-header";
    header.innerHTML = `
      <div class="object-title-container">
        <div class="object-icon">ℹ️</div>
        <div class="object-name">${message}</div>
      </div>
      <button class="object-close" onclick="this.parentElement.parentElement.remove(); window.updateScrollHandles && window.updateScrollHandles();">
        <img src="close.svg" alt="Close" width="16" height="16" />
      </button>
    `;

    column.appendChild(header);
    gldiv.appendChild(column);

    if (window.updateScrollHandles) {
      window.updateScrollHandles();
    }

    return column;
  }

  // Create loading column for analysis
  createLoadingColumn(message, croppedImageBase64 = null) {
    const gldiv = document.getElementById("gldiv");
    if (!gldiv) return null;

    const loadingColumn = document.createElement("div");
    loadingColumn.className = "object-column";

    // Create header
    const header = document.createElement("div");
    header.className = "object-header";
    header.innerHTML = `
      <div class="object-title-container">
        <div class="object-icon">🔬</div>
        <div class="object-name">Analyzing...</div>
      </div>
    `;

    loadingColumn.appendChild(header);

    // Add cropped image if provided
    if (croppedImageBase64) {
      const imageContainer = document.createElement("div");
      imageContainer.className = "cropped-image-container";
      
      const img = document.createElement("img");
      img.src = `data:image/jpeg;base64,${croppedImageBase64}`;
      img.className = "image-highlighted";
      img.alt = "Analysis region";
      
      imageContainer.appendChild(img);
      loadingColumn.appendChild(imageContainer);
    }

    // Add loading indicator
    const loadingIndicator = document.createElement("div");
    loadingIndicator.className = "loading-indicator loading-message-italic";
    loadingIndicator.textContent = "Processing...";
    loadingColumn.appendChild(loadingIndicator);

    gldiv.appendChild(loadingColumn);

    // Update scroll handles if available
    if (window.updateScrollHandles) {
      window.updateScrollHandles();
    }

    return loadingColumn;
  }

  // Show main app interface
  showMainApp() {
    const mainInterface = document.getElementById('main-app-interface');
    if (mainInterface) {
      mainInterface.style.display = 'block';
    }
  }

  // Remove development blur/payment states (for HTTPS) - but preserve valid payment modal
  clearDevelopmentStates() {
    if (location.protocol === 'https:') {
      
      
      setTimeout(() => {
        const mainInterface = document.getElementById('main-app-interface');
        const paymentModal = document.getElementById('payment-modal');
        
        if (mainInterface) {
          // Only clear payment-required if modal is not legitimately shown
          const modalIsShown = paymentModal && paymentModal.classList.contains('show');
          
          if (!modalIsShown) {
            mainInterface.classList.remove('payment-required');
          }
          
          mainInterface.style.filter = 'none';
          mainInterface.style.opacity = '1';
          mainInterface.style.pointerEvents = 'auto';
          mainInterface.style.visibility = 'visible';
          mainInterface.style.display = 'block';
        }
        
        ');
      }, 100); // Reduced timeout to minimize flash
    }
  }

  // File conversion utilities
  fileToBase64(file) {
    return new Promise((resolve, reject) => {
      const reader = new FileReader();
      reader.onload = () => {
        const result = reader.result.split(",")[1];
        resolve(result);
      };
      reader.onerror = reject;
      reader.readAsDataURL(file);
    });
  }

  async urlToBase64(url) {
    try {
      const response = await fetch(url);
      const blob = await response.blob();
      
      return new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = () => {
          const result = reader.result.split(",")[1];
          resolve(result);
        };
        reader.onerror = reject;
        reader.readAsDataURL(blob);
      });
    } catch (error) {
      throw new Error(`Failed to fetch image from URL: ${error.message}`);
    }
  }

  // Molecule name helper
  getMoleculeName(chemical) {
    const moleculeNames = {
      O: "Water",
      CCO: "Ethanol", 
      "CC(=O)O": "Acetic Acid",
      C: "Methane",
      CO: "Methanol",
      "C(CO)N": "Ethanolamine",
      "C(C(=O)O)N": "Glycine",
      "C(CC(=O)O)N": "GABA",
      "C(CC(=O)O)C(=O)O": "Succinic Acid",
      "C(C(=O)O)O": "Glycolic Acid",
      "C1=CC=CC=C1": "Benzene",
      "CC1=CC=CC=C1": "Toluene",
      "C1=CC=C(C=C1)O": "Phenol",
      "C1=CC=C(C=C1)N": "Aniline",
      "CN1C=NC2=C1C(=O)N(C(=O)N2C)C": "Caffeine",
      "CC(=O)OC1=CC=CC=C1C(=O)O": "Aspirin"
    };

    return moleculeNames[chemical.smiles] || chemical.name || "Unknown Molecule";
  }

  // Format result message
  createResultMessage(icon, objectName, smilesCount, useQuotes = false) {
    const name = useQuotes ? `"${objectName}"` : objectName;
    const plural = smilesCount !== 1 ? "s" : "";
    return `${icon} ${name} → ${smilesCount} molecule${plural} found`;
  }

  // Development debugging functions
  setupDebuggingFunctions() {
    // Payment bypass toggle
    window.bypassPayment = () => {
      const mainInterface = document.getElementById('main-app-interface');
      const paymentPopdown = document.getElementById('payment-modal');
      
      if (mainInterface.classList.contains('payment-required')) {
        
        mainInterface.classList.remove('payment-required');
        mainInterface.style.opacity = '1';
        mainInterface.style.filter = 'none';
        mainInterface.style.pointerEvents = 'auto';
        if (paymentPopdown) {
          paymentPopdown.style.display = 'none';
        }
      } else {
        
        mainInterface.classList.add('payment-required');
        mainInterface.style.opacity = '';
        mainInterface.style.filter = '';
        mainInterface.style.pointerEvents = '';
        if (paymentPopdown) {
          paymentPopdown.style.display = 'block';
        }
      }
    };

    // Enable developer mode
    window.enableDevMode = () => {
      
      
      // Set up developer account
      if (window.paymentManager) {
        paymentManager.setupDeveloperAccount();
        
        
        // Update app payment setup status
        if (window.app) {
          app.hasPaymentSetup = true;
          
        }
        
        // Hide payment modal
        const paymentPopdown = document.getElementById('payment-modal');
        if (paymentPopdown) {
          paymentManager.hidePaymentModal();
        }
        
        
      } else {
        console.error('❌ Payment manager not available');
      }
    };

    // Show both sections clearly
    window.showBothSections = () => {
      
      const mainInterface = document.getElementById('main-app-interface');
      const paymentPopdown = document.getElementById('payment-modal');
      
      if (paymentPopdown) {
        paymentPopdown.style.display = 'block';
        paymentPopdown.style.visibility = 'visible';
      }
      
      if (mainInterface) {
        mainInterface.style.display = 'block';
        mainInterface.style.visibility = 'visible';
        mainInterface.style.opacity = '1';
        mainInterface.style.filter = 'none';
        mainInterface.style.pointerEvents = 'auto';
        mainInterface.classList.remove('payment-required');
      }
      
      
    };

    // Simple app show
    window.showApp = () => {
      
      const mainInterface = document.getElementById('main-app-interface');
      
      if (mainInterface) {
        mainInterface.style.display = 'block';
        mainInterface.style.opacity = '1';
        mainInterface.style.filter = 'none';
        mainInterface.style.pointerEvents = 'auto';
        mainInterface.style.visibility = 'visible';
        mainInterface.classList.remove('payment-required');
        
      } else {
        console.error('❌ Main interface element not found');
      }
    };

    // Element checker
    window.checkElements = () => {
      
      );
      );
      );
      );
    };
  }

  // Cleanup function
  cleanup() {
    // Remove event listeners and clean up resources
    if (this.cameraMode) {
      this.cameraMode.removeEventListener("change", this.updateInputMode);
    }
    if (this.photoMode) {
      this.photoMode.removeEventListener("change", this.updateInputMode);
    }
  }
}

// Export singleton instance
export const uiManager = new UIManager(); 