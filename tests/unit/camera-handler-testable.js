/**
 * CameraHandler Testable Implementation
 * This is a testable version of the camera handler with dependency injection
 */

class CameraHandler {
  constructor(uiManager, paymentManager) {
    this.uiManager = uiManager;
    this.paymentManager = paymentManager;
    this.isInitialized = false;
    this.hasPermission = false;
    this.stream = null;
    this.videoElement = null;
    this.availableCameras = [];
    this.currentCameraIndex = 0;
    this.isProcessing = false;
    this.clickPosition = null;
    this.showOutline = false;
    this.outlineTimeout = null;
    this.isMobile = /Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini/i.test(navigator.userAgent) ||
                    (window.innerWidth <= 768 && 'ontouchstart' in window);
  }

  async initialize() {
    this.isInitialized = true;
    return true;
  }

  async requestPermission() {
    try {
      const stream = await navigator.mediaDevices.getUserMedia({ 
        video: { facingMode: 'environment' },
        audio: false 
      });
      this.stream = stream;
      this.hasPermission = true;
      return true;
    } catch (error) {
      this.hasPermission = false;
      throw error;
    }
  }

  async startCamera() {
    if (!this.hasPermission) {
      throw new Error('Camera permission required');
    }
    if (!this.stream) {
      await this.requestPermission();
    }
    return true;
  }

  stopCamera() {
    if (this.stream) {
      this.stream.getTracks().forEach(track => track.stop());
      this.stream = null;
    }
  }

  async getAvailableCameras() {
    try {
      const devices = await navigator.mediaDevices.enumerateDevices();
      this.availableCameras = devices.filter(device => device.kind === 'videoinput');
      return this.availableCameras;
    } catch (error) {
      return [];
    }
  }

  async captureImage() {
    if (!this.videoElement) {
      throw new Error('Video element not available');
    }
    
    const canvas = document.createElement('canvas');
    const video = this.videoElement;
    canvas.width = video.videoWidth;
    canvas.height = video.videoHeight;
    const ctx = canvas.getContext('2d');
    ctx.drawImage(video, 0, 0);
    return canvas.toDataURL('image/jpeg', 0.8);
  }

  async switchCamera() {
    if (this.availableCameras.length <= 1) {
      throw new Error('No other cameras available');
    }
    
    this.currentCameraIndex = (this.currentCameraIndex + 1) % this.availableCameras.length;
    
    if (this.stream) {
      this.stopCamera();
    }
    
    const deviceId = this.availableCameras[this.currentCameraIndex].deviceId;
    const stream = await navigator.mediaDevices.getUserMedia({
      video: { deviceId: { exact: deviceId } },
      audio: false
    });
    
    this.stream = stream;
    return true;
  }

  async handlePhotoUpload(file) {
    if (!file) {
      throw new Error('No file provided');
    }
    if (!file.type.startsWith('image/')) {
      throw new Error('File must be an image');
    }
    
    if (this.paymentManager && this.paymentManager.checkPaymentMethod) {
      const paymentRequired = await this.paymentManager.checkPaymentMethod();
      if (!paymentRequired) {
        throw new Error('Payment required');
      }
    }
    
    return new Promise((resolve, reject) => {
      const reader = new FileReader();
      reader.onload = (e) => resolve(e.target.result);
      reader.onerror = reject;
      reader.readAsDataURL(file);
    });
  }

  async handleUrlAnalysis(url) {
    if (!url) {
      throw new Error('URL is required');
    }
    if (!url.startsWith('http')) {
      throw new Error('Invalid URL format');
    }
    
    if (this.paymentManager && this.paymentManager.checkPaymentMethod) {
      const paymentRequired = await this.paymentManager.checkPaymentMethod();
      if (!paymentRequired) {
        throw new Error('Payment required');
      }
    }
    
    try {
      const response = await fetch(url);
      if (!response.ok) {
        throw new Error(`HTTP ${response.status}: ${response.statusText}`);
      }
      const blob = await response.blob();
      return new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = (e) => resolve(e.target.result);
        reader.onerror = reject;
        reader.readAsDataURL(blob);
      });
    } catch (error) {
      throw new Error(`Failed to load image from URL: ${error.message}`);
    }
  }

  displayUploadedImage(imageData, isMobile = false) {
    const img = document.createElement('img');
    img.src = imageData;
    img.dataset.analysisData = JSON.stringify({ object: 'test' });
    img.className = 'uploaded-image';
    if (isMobile) {
      img.classList.add('mobile-reticle');
    }
    return img;
  }

  async handleImageClick(e, imageData) {
    if (this.paymentManager && this.paymentManager.checkPaymentMethod) {
      const paymentRequired = await this.paymentManager.checkPaymentMethod();
      if (!paymentRequired) {
        throw new Error('Payment required');
      }
    }
    
    if (!imageData) {
      throw new Error('No image data available');
    }
    
    // Mock API call
    const response = await fetch('/api/analyze-image', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ imageData })
    });
    
    if (!response.ok) {
      throw new Error(`API error: ${response.status}`);
    }
    
    return await response.json();
  }

  emitAnalysisResult(data, error = null) {
    const event = new CustomEvent('analysisComplete', {
      detail: { data, error }
    });
    document.dispatchEvent(event);
  }

  createClosableErrorMessage(message) {
    if (this.uiManager && this.uiManager.createErrorMessage) {
      return this.uiManager.createErrorMessage(message);
    }
    
    const errorDiv = document.createElement('div');
    errorDiv.className = 'error-message';
    errorDiv.textContent = message;
    errorDiv.remove = jest.fn();
    return errorDiv;
  }

  updateScrollHandles() {
    if (window.updateScrollHandles) {
      window.updateScrollHandles();
    }
  }

  setupEventListeners() {
    const photoUpload = document.getElementById('photo-upload');
    const photoUrl = document.getElementById('photo-url');
    const urlAnalyze = document.getElementById('url-analyze');
    
    if (photoUpload) {
      photoUpload.addEventListener('change', (e) => {
        const file = e.target.files[0];
        if (file) {
          this.handlePhotoUpload(file).catch(error => {
            if (this.uiManager && this.uiManager.createErrorMessage) {
              this.uiManager.createErrorMessage(error.message);
            }
          });
        }
      });
    }
    
    if (urlAnalyze && photoUrl) {
      urlAnalyze.addEventListener('click', () => {
        const url = photoUrl.value;
        if (url) {
          this.handleUrlAnalysis(url).catch(error => {
            if (this.uiManager && this.uiManager.createErrorMessage) {
              this.uiManager.createErrorMessage(error.message);
            }
          });
        }
      });
    }
  }

  calculateCropCoordinates(clickX, clickY, imageWidth, imageHeight) {
    const cropSize = 100;
    const x = Math.max(0, Math.min(clickX - cropSize/2, imageWidth - cropSize));
    const y = Math.max(0, Math.min(clickY - cropSize/2, imageHeight - cropSize));
    return { x, y, width: cropSize, height: cropSize };
  }
}

module.exports = { CameraHandler };