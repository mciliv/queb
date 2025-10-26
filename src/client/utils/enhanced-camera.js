// Enhanced Camera APIs for PWA
// Provides better camera access with fallbacks and advanced features
import logger from '../core/logger.js';

class EnhancedCamera {
  constructor() {
    this.stream = null;
    this.videoTrack = null;
    this.capabilities = null;
    this.settings = null;
  }

  /**
   * Get enhanced camera capabilities
   */
  async getCapabilities() {
    if (!this.videoTrack) {
      throw new Error('No active camera stream');
    }

    this.capabilities = this.videoTrack.getCapabilities();
    this.settings = this.videoTrack.getSettings();
    
    return {
      capabilities: this.capabilities,
      settings: this.settings,
      supportedConstraints: navigator.mediaDevices.getSupportedConstraints(),
      hasAdvancedFeatures: !!(this.capabilities.focusMode || this.capabilities.exposureMode)
    };
  }

  /**
   * Request camera with enhanced constraints
   */
  async requestCamera(options = {}) {
    const defaultConstraints = {
      video: {
        facingMode: options.facingMode || 'environment',
        width: { ideal: 1920, min: 720 },
        height: { ideal: 1080, min: 480 },
        aspectRatio: { ideal: 16/9 },
        frameRate: { ideal: 30, min: 15 }
      }
    };

    // Add advanced constraints if supported
    if (navigator.mediaDevices.getSupportedConstraints().focusMode) {
      defaultConstraints.video.focusMode = 'continuous';
    }

    if (navigator.mediaDevices.getSupportedConstraints().exposureMode) {
      defaultConstraints.video.exposureMode = 'continuous';
    }

    if (navigator.mediaDevices.getSupportedConstraints().whiteBalanceMode) {
      defaultConstraints.video.whiteBalanceMode = 'continuous';
    }

    try {
      this.stream = await navigator.mediaDevices.getUserMedia(defaultConstraints);
      this.videoTrack = this.stream.getVideoTracks()[0];
      
      // Get capabilities after successful stream
      const capabilities = await this.getCapabilities();
      
      return {
        stream: this.stream,
        videoTrack: this.videoTrack,
        ...capabilities
      };
    } catch (error) {
      logger.warn('Enhanced camera request failed, using fallback', { error: error.message });

      // Fallback to basic constraints
      const fallbackConstraints = {
        video: {
          facingMode: options.facingMode || 'environment'
        }
      };
      
      this.stream = await navigator.mediaDevices.getUserMedia(fallbackConstraints);
      this.videoTrack = this.stream.getVideoTracks()[0];
      
      return {
        stream: this.stream,
        videoTrack: this.videoTrack,
        capabilities: this.videoTrack.getCapabilities(),
        settings: this.videoTrack.getSettings(),
        supportedConstraints: navigator.mediaDevices.getSupportedConstraints(),
        hasAdvancedFeatures: false
      };
    }
  }

  /**
   * Switch camera (front/back)
   */
  async switchCamera() {
    if (!this.videoTrack) {
      throw new Error('No active camera stream');
    }

    const currentFacingMode = this.settings?.facingMode;
    const newFacingMode = currentFacingMode === 'user' ? 'environment' : 'user';
    
    // Stop current stream
    this.stopCamera();
    
    // Request new camera
    return await this.requestCamera({ facingMode: newFacingMode });
  }

  /**
   * Apply camera settings (focus, exposure, etc.)
   */
  async applySettings(settings = {}) {
    if (!this.videoTrack) {
      throw new Error('No active camera stream');
    }

    const constraints = {};
    
    // Focus mode
    if (settings.focusMode && this.capabilities?.focusMode?.includes(settings.focusMode)) {
      constraints.focusMode = settings.focusMode;
    }
    
    // Exposure mode  
    if (settings.exposureMode && this.capabilities?.exposureMode?.includes(settings.exposureMode)) {
      constraints.exposureMode = settings.exposureMode;
    }
    
    // White balance
    if (settings.whiteBalanceMode && this.capabilities?.whiteBalanceMode?.includes(settings.whiteBalanceMode)) {
      constraints.whiteBalanceMode = settings.whiteBalanceMode;
    }

    // Torch (flashlight)
    if (settings.torch !== undefined && this.capabilities?.torch) {
      constraints.torch = settings.torch;
    }

    if (Object.keys(constraints).length > 0) {
      try {
        await this.videoTrack.applyConstraints({ advanced: [constraints] });
        this.settings = this.videoTrack.getSettings();
        return this.settings;
      } catch (error) {
        logger.warn('Failed to apply camera settings', { error: error.message });
        return this.settings;
      }
    }
    
    return this.settings;
  }

  /**
   * Capture high-quality image
   */
  async captureImage(videoElement, options = {}) {
    if (!videoElement || !this.videoTrack) {
      throw new Error('No video element or camera stream');
    }

    const canvas = document.createElement('canvas');
    const ctx = canvas.getContext('2d');
    
    // Use video dimensions or specified size
    canvas.width = options.width || videoElement.videoWidth || 1920;
    canvas.height = options.height || videoElement.videoHeight || 1080;
    
    // Draw video frame to canvas
    ctx.drawImage(videoElement, 0, 0, canvas.width, canvas.height);
    
    // Convert to blob with high quality
    return new Promise((resolve) => {
      canvas.toBlob(resolve, options.format || 'image/jpeg', options.quality || 0.95);
    });
  }

  /**
   * Get available cameras
   */
  async getAvailableCameras() {
    const devices = await navigator.mediaDevices.enumerateDevices();
    return devices
      .filter(device => device.kind === 'videoinput')
      .map(device => ({
        deviceId: device.deviceId,
        label: device.label || `Camera ${device.deviceId.slice(0, 8)}`,
        facingMode: this.guessFacingMode(device.label)
      }));
  }

  /**
   * Guess facing mode from device label
   */
  guessFacingMode(label) {
    const lowerLabel = label.toLowerCase();
    if (lowerLabel.includes('front') || lowerLabel.includes('user')) {
      return 'user';
    }
    if (lowerLabel.includes('back') || lowerLabel.includes('rear') || lowerLabel.includes('environment')) {
      return 'environment';
    }
    return 'unknown';
  }

  /**
   * Stop camera stream
   */
  stopCamera() {
    if (this.stream) {
      this.stream.getTracks().forEach(track => track.stop());
      this.stream = null;
      this.videoTrack = null;
      this.capabilities = null;
      this.settings = null;
    }
  }

  /**
   * Check if device supports advanced camera features
   */
  static hasAdvancedCameraSupport() {
    const constraints = navigator.mediaDevices?.getSupportedConstraints() || {};
    return !!(constraints.focusMode || constraints.exposureMode || constraints.torch);
  }

  /**
   * Check if running as installed PWA
   */
  static isPWA() {
    return window.matchMedia('(display-mode: standalone)').matches ||
           window.navigator.standalone === true;
  }
}

export default EnhancedCamera;
