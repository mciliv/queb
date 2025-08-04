// InputSection.js - Handles text input, camera, and photo upload modes
export class InputSection {
  constructor() {
    this.objectInput = null;
    this.cameraMode = null;
    this.photoMode = null;
    this.cameraContainer = null;
    this.photoOptions = null;
  }

  initialize() {
    this.objectInput = document.getElementById("object-input");
    this.cameraMode = document.getElementById("camera-mode");
    this.photoMode = document.getElementById("photo-mode");
    this.cameraContainer = document.getElementById("camera-container");
    this.photoOptions = document.getElementById("photo-options");

    this.setupEventListeners();
    console.log('✅ Input section initialized');
  }

  setupEventListeners() {
    // Mode switching
    if (this.cameraMode) {
      this.cameraMode.addEventListener('change', () => this.handleModeChange());
    }

    if (this.photoMode) {
      this.photoMode.addEventListener('change', () => this.handleModeChange());
    }

    // Photo upload handling
    const photoUpload = document.getElementById('photo-upload');
    if (photoUpload) {
      photoUpload.addEventListener('change', (e) => this.handlePhotoUpload(e));
    }

    // URL analysis
    const urlAnalyze = document.getElementById('url-analyze');
    if (urlAnalyze) {
      urlAnalyze.addEventListener('click', () => this.handleUrlAnalysis());
    }

    // URL input enter key
    const photoUrl = document.getElementById('photo-url');
    if (photoUrl) {
      photoUrl.addEventListener('keyup', (e) => {
        if (e.key === 'Enter') {
          this.handleUrlAnalysis();
        }
      });
    }
  }

  handleModeChange() {
    const isCameraMode = this.cameraMode?.checked;
    const isPhotoMode = this.photoMode?.checked;

    // Show/hide camera container
    if (this.cameraContainer) {
      this.cameraContainer.style.display = isCameraMode ? 'block' : 'none';
    }

    // Show/hide photo options
    if (this.photoOptions) {
      this.photoOptions.style.display = isPhotoMode ? 'block' : 'none';
    }

    // Dispatch mode change event
    document.dispatchEvent(new CustomEvent('inputModeChanged', {
      detail: { cameraMode: isCameraMode, photoMode: isPhotoMode }
    }));
  }

  async handlePhotoUpload(event) {
    const file = event.target.files[0];
    if (!file) return;

    try {
      const imageData = await this.readFileAsDataURL(file);
      await this.analyzeImage(imageData, file.name);
    } catch (error) {
      console.error('❌ Photo upload failed:', error);
      document.dispatchEvent(new CustomEvent('appError', {
        detail: { message: 'Failed to upload photo', type: 'error' }
      }));
    }
  }

  async handleUrlAnalysis() {
    const photoUrl = document.getElementById('photo-url');
    if (!photoUrl || !photoUrl.value.trim()) return;

    try {
      await this.analyzeImage(photoUrl.value.trim(), 'URL Image');
      photoUrl.value = ''; // Clear input after analysis
    } catch (error) {
      console.error('❌ URL analysis failed:', error);
      document.dispatchEvent(new CustomEvent('appError', {
        detail: { message: 'Failed to analyze image URL', type: 'error' }
      }));
    }
  }

  readFileAsDataURL(file) {
    return new Promise((resolve, reject) => {
      const reader = new FileReader();
      reader.onload = (e) => resolve(e.target.result);
      reader.onerror = (e) => reject(e);
      reader.readAsDataURL(file);
    });
  }

  async analyzeImage(imageData, imageName) {
    try {
      const response = await fetch('/api/analyze-image', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ 
          image: imageData,
          imageName: imageName
        })
      });

      if (!response.ok) throw new Error('Image analysis failed');
      
      const result = await response.json();
      
      // Dispatch analysis complete event
      document.dispatchEvent(new CustomEvent('imageAnalysisComplete', {
        detail: {
          output: result.output,
          icon: null,
          objectName: imageName,
          useQuotes: false,
          croppedImageData: imageData
        }
      }));
    } catch (error) {
      throw new Error(`Image analysis failed: ${error.message}`);
    }
  }

  // Public methods
  focusInput() {
    if (this.objectInput) {
      this.objectInput.focus();
      this.objectInput.select();
    }
  }

  clearInput() {
    if (this.objectInput) {
      this.objectInput.value = '';
    }
  }

  setInputValue(value) {
    if (this.objectInput) {
      this.objectInput.value = value;
    }
  }

  getInputValue() {
    return this.objectInput ? this.objectInput.value.trim() : '';
  }

  enableInput() {
    if (this.objectInput) {
      this.objectInput.disabled = false;
      this.objectInput.style.opacity = '1';
    }
  }

  disableInput() {
    if (this.objectInput) {
      this.objectInput.disabled = true;
      this.objectInput.style.opacity = '0.5';
    }
  }

  // Mode management
  setCameraMode(enabled) {
    if (this.cameraMode) {
      this.cameraMode.checked = enabled;
      this.handleModeChange();
    }
  }

  setPhotoMode(enabled) {
    if (this.photoMode) {
      this.photoMode.checked = enabled;
      this.handleModeChange();
    }
  }

  getCurrentMode() {
    return {
      camera: this.cameraMode?.checked || false,
      photo: this.photoMode?.checked || false
    };
  }
} 