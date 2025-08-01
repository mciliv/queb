import { paymentManager } from './payment.js';
import { uiManager } from './ui-utils.js';
import { logger } from './logger.js';

class CameraManager {
  constructor() {
    this.video = null;
    this.currentStream = null;
    this.facingMode = "environment";
    
    this.permissionKey = 'mol_camera_permission';
    this.deviceId = this.getDeviceId();
    this.isMobile = /Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini/i.test(navigator.userAgent);
    this.isSafari = /^((?!chrome|android).)*safari/i.test(navigator.userAgent);
    this.isIOS = /iPad|iPhone|iPod/.test(navigator.userAgent);
    this.cameraPermission = this.hasStoredCameraPermission();
  }

  getDeviceId() {
    let deviceId = localStorage.getItem('mol_device_id');
    if (!deviceId) {
      deviceId = 'device_' + Date.now() + '_' + Math.random().toString(36).substr(2, 9);
      localStorage.setItem('mol_device_id', deviceId);
    }
    return deviceId;
  }

  saveCameraPermission(granted) {
    const permissions = this.getStoredPermissions();
    permissions[this.deviceId] = {
      camera: granted,
      timestamp: Date.now()
    };
    localStorage.setItem(this.permissionKey, JSON.stringify(permissions));
  }

  // Get stored permissions
  getStoredPermissions() {
    try {
      const stored = localStorage.getItem(this.permissionKey);
      return stored ? JSON.parse(stored) : {};
    } catch (error) {
      return {};
    }
  }

  // Check if camera permission was previously granted
  hasStoredCameraPermission() {
    const permissions = this.getStoredPermissions();
    const devicePerms = permissions[this.deviceId];
    if (!devicePerms) return false;
    
    // Check if permission is still valid (not older than 30 days)
    const thirtyDaysAgo = Date.now() - (30 * 24 * 60 * 60 * 1000);
    return devicePerms.camera && devicePerms.timestamp > thirtyDaysAgo;
  }

  // Clear stored camera permission
  clearStoredPermission() {
    const permissions = this.getStoredPermissions();
    delete permissions[this.deviceId];
    localStorage.setItem(this.permissionKey, JSON.stringify(permissions));
  }

  // Get permission status for UI display
  getPermissionStatus() {
    if (this.hasStoredCameraPermission()) {
      return 'granted';
    } else if (this.currentStream) {
      return 'active';
    } else {
      return 'denied';
    }
  }

  // Check if camera is ready for interaction
  isCameraReady() {
    return this.video && 
           this.currentStream && 
           this.video.videoWidth > 0 && 
           this.video.videoHeight > 0 && 
           this.video.readyState >= 2;
  }

  // Update UI based on permission status
  updatePermissionUI() {
    const status = this.getPermissionStatus();
    const video = this.video;
    
    if (status === 'granted' || status === 'active') {
      video.classList.add('permission-granted');
      video.classList.remove('permission-denied');
    } else {
      video.classList.add('permission-denied');
      video.classList.remove('permission-granted');
    }
  }

  // Show permission request button
  showPermissionRequest() {
    const video = this.video;
    const permissionBtn = document.createElement('button');
    permissionBtn.className = 'camera-permission-btn';
    permissionBtn.textContent = 'Enable Camera';
    permissionBtn.onclick = () => this.requestPermission();
    
    video.parentNode.appendChild(permissionBtn);
  }

  // Initialize camera system
  async initialize() {
    this.video = document.getElementById("video-feed");
    
    if (!navigator.mediaDevices || !navigator.mediaDevices.getUserMedia) {
      this.showError("Camera API not supported in this browser");
      return false;
    }

    const isSecureContext = window.isSecureContext || location.protocol === "https:";
    if (!isSecureContext && !location.hostname.includes("localhost")) {
      this.showError("Camera requires HTTPS on mobile devices. Please use HTTPS or localhost.");
      return false;
    }

    this.setupEventListeners();
    this.updateCameraInstructions();
    
    // Check if we have stored permission and auto-start camera
    if (this.hasStoredCameraPermission()) {
      try {
        await this.startCamera();
        await this.setupSwitchCamera();
        this.updateReticleVisibility();
        return true;
      } catch (err) {
        // If auto-start fails, clear the stored permission
        this.saveCameraPermission(false);
        this.updateCameraInstructions();
      }
    }
    
    return false;
  }

  // Setup camera event listeners
  setupEventListeners() {
    this.video.addEventListener("click", (e) => this.handleInteraction(e));
    this.video.addEventListener("touchstart", (e) => {
      e.preventDefault();
      this.handleInteraction(e.touches[0]);
    });

    // Safari-specific touch handling
    if (this.isSafari || this.isIOS) {
      this.video.addEventListener("touchend", (e) => e.preventDefault());
      
      let lastTouchEnd = 0;
      this.video.addEventListener("touchend", (e) => {
        const now = new Date().getTime();
        if (now - lastTouchEnd <= 300) {
          e.preventDefault();
        }
        lastTouchEnd = now;
      }, false);
    }
  }

  // Start camera stream
  async startCamera() {
    if (paymentManager.isPaymentRequired()) {
      return;
    }

    this.stopCurrentStream();
    
    // Validate video element
    if (!this.video) {
      logger.error('Video element not found');
      return;
    }

    try {
      let stream;

      if (this.isSafari || this.isIOS) {
        try {
          stream = await navigator.mediaDevices.getUserMedia(this.getSafariConstraints());
        } catch (err) {
          stream = await navigator.mediaDevices.getUserMedia(this.getBasicConstraints());
        }
      } else {
        try {
          stream = await navigator.mediaDevices.getUserMedia(this.getSimpleConstraints());
        } catch (err) {
          stream = await navigator.mediaDevices.getUserMedia(this.getBasicConstraints());
        }
      }

      this.currentStream = stream;
      this.video.srcObject = stream;

      // Safari-specific video attributes
      this.video.setAttribute("playsinline", "true");
      this.video.setAttribute("webkit-playsinline", "true");
      this.video.setAttribute("x-webkit-airplay", "allow");

      await this.video.play();
      this.hideError();
      this.updateReticleVisibility();
      
      // Save successful camera permission
      this.saveCameraPermission(true);
      
      logger.info('Camera started successfully', {
        videoWidth: this.video.videoWidth,
        videoHeight: this.video.videoHeight,
        readyState: this.video.readyState,
        srcObject: !!this.video.srcObject
      });
      
      // Update instruction text based on camera readiness
      this.updateCameraInstructions();
      
    } catch (err) {
      this.saveCameraPermission(false);
      this.handleCameraError(err);
    }
  }

  // Stop current camera stream
  stopCurrentStream() {
    if (this.currentStream) {
      this.currentStream.getTracks().forEach(track => track.stop());
      this.currentStream = null;
    }
  }

  // Get camera constraints for different browsers
  getSimpleConstraints() {
    return {
      video: {
        facingMode: this.facingMode,
        width: { ideal: 1280 },
        height: { ideal: 720 },
      },
    };
  }

  getBasicConstraints() {
    return {
      video: {
        width: { ideal: 1280 },
        height: { ideal: 720 },
      },
    };
  }

  getSafariConstraints() {
    return {
      video: {
        facingMode: this.facingMode,
        width: { min: 640, ideal: 1280, max: 1920 },
        height: { min: 480, ideal: 720, max: 1080 },
      },
    };
  }

  // Handle camera interaction (click/tap)
  async handleInteraction(evt) {
    logger.cameraEvent('interaction_detected', { 
      isMobile: this.isMobile,
      eventType: evt.type 
    });
    
    let paymentCheck = false;
    try {
      paymentCheck = await this.checkPaymentSetupForSidebar();
      logger.debug('Payment check result', { paymentCheck });
    } catch (error) {
      logger.warn('Payment check failed, proceeding anyway', error);
      paymentCheck = true; // Fallback to allow analysis
    }
    
    if (!paymentCheck) {
      logger.warn('Payment required - showing message');
      // Show simple message directing to sidebar
      const messageColumn = uiManager.createColumn("Payment setup required", "payment-required");
      messageColumn.innerHTML = `
        <div class="molecule-container">
          <div class="molecule-info">
            <h3>Payment Setup Required</h3>
            <p>Complete payment setup in the sidebar to analyze molecules via camera</p>
          </div>
        </div>
      `;
      return;
    }

    logger.info('Payment check passed, proceeding with analysis');

    // Mobile reticle validation
    if (this.isMobile && !this.isWithinReticle(evt)) {
      logger.warn('Mobile reticle validation failed');
      this.showReticleFeedback();
      return;
    }

    logger.info('Proceeding with camera analysis');
    
    // Check if camera is ready
    if (!this.isCameraReady()) {
      logger.warn('Camera not ready, attempting to restart');
      try {
        await this.startCamera();
        // Wait a moment for camera to stabilize
        await new Promise(resolve => setTimeout(resolve, 1000));
        
        if (!this.isCameraReady()) {
          throw new Error('Camera failed to initialize - please refresh the page and try again');
        }
      } catch (restartError) {
        throw new Error('Camera initialization failed - please check camera permissions and try again');
      }
    }
    
    this.showCropOutline(evt);
    
    try {
      // Validate camera state before capture
      if (!this.video || !this.currentStream) {
        throw new Error('Camera not ready - please wait for camera to initialize');
      }
      
      // Capture and analyze the image
      logger.cameraEvent('capturing_image');
      const captureData = await this.captureAndAnalyze(evt);
      logger.info('Image captured successfully');
      
      // Create loading column
      const loadingColumn = uiManager.createLoadingColumn("Analyzing...", captureData.croppedBase64);
      
      // Send to server for analysis
      logger.analysisEvent('camera_analysis_request', {
        coordinates: captureData.coordinates,
        hasImageData: !!captureData.imageBase64
      });
      
      // Create AbortController for timeout
      const controller = new AbortController();
      const timeoutId = setTimeout(() => controller.abort(), 45000); // 45 second timeout
      
      const response = await fetch("/image-molecules", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          imageBase64: captureData.imageBase64,
          croppedImageBase64: captureData.croppedBase64,
          x: captureData.coordinates.x,
          y: captureData.coordinates.y,
          cropMiddleX: captureData.coordinates.cropMiddleX,
          cropMiddleY: captureData.coordinates.cropMiddleY,
          cropSize: captureData.coordinates.cropSize,
        }),
        signal: controller.signal
      });

      clearTimeout(timeoutId);

      logger.debug('Camera analysis response status', { status: response.status });
      if (!response.ok) {
        const errorText = await response.text();
        throw new Error(`HTTP ${response.status}: ${errorText}`);
      }
      
      const { output } = await response.json();
      logger.analysisEvent('camera_analysis_completed', { output });

      // Remove loading column
      loadingColumn.remove();
      
      // Emit analysis result event for app to handle
      const event = new CustomEvent('imageAnalysisComplete', {
        detail: { 
          output, 
          icon: "📷", 
          objectName: output.object || "Camera capture", 
          useQuotes: false, 
          croppedImageData: captureData.croppedBase64 
        }
      });
      document.dispatchEvent(event);
      
      // Try to increment usage, but don't fail if it doesn't work
      try {
        await this.incrementUsageForSidebar();
      } catch (usageError) {
        logger.warn('Usage increment failed', usageError);
      }
      
    } catch (error) {
      logger.error('Camera analysis error', error);
      
      // Enhanced error handling
      let errorMessage = 'Analysis failed: ';
      
      if (error.name === 'AbortError') {
        errorMessage = 'Camera analysis timed out. Please check your internet connection and try again.';
      } else if (error.message.includes('Failed to fetch') || error.message.includes('fetch')) {
        errorMessage = 'Network connection failed. Please check your internet connection and try again.';
      } else if (error.message.includes('Network connection failed')) {
        errorMessage = 'Unable to connect to analysis service. Please check your internet connection.';
      } else if (error.message.includes('Camera not ready')) {
        errorMessage = error.message; // Keep camera-specific messages as-is
      } else {
        errorMessage += error.message;
      }
      
      this.createClosableErrorMessage(errorMessage);
    }
  }

  // Check if tap is within mobile reticle
  isWithinReticle(evt) {
    const rect = this.video.getBoundingClientRect();
    const centerX = rect.left + rect.width / 2;
    const centerY = rect.top + rect.height / 2;
    const tapX = evt.clientX;
    const tapY = evt.clientY;

    const reticleRadius = 50;
    const distanceFromCenter = Math.sqrt(
      Math.pow(tapX - centerX, 2) + Math.pow(tapY - centerY, 2)
    );

    return distanceFromCenter <= reticleRadius;
  }

  // Show reticle feedback animation
  showReticleFeedback() {
    const reticle = document.querySelector(".mobile-reticle");
    if (reticle) {
      reticle.style.animation = "reticlePulse 0.5s ease";
      setTimeout(() => {
        reticle.style.animation = "";
      }, 500);
    }
  }

  // Show crop outline at interaction point
  showCropOutline(evt) {
    const cropSize = 100;
    const outline = document.createElement("div");
    outline.className = "crop-outline";
    outline.style.width = cropSize + "px";
    outline.style.height = cropSize + "px";
    
    // Use page coordinates for fixed positioning
    const x = evt.pageX || evt.clientX;
    const y = evt.pageY || evt.clientY;
    
    outline.style.left = (x - cropSize / 2) + "px";
    outline.style.top = (y - cropSize / 2) + "px";
    
    document.body.appendChild(outline);
    
    setTimeout(() => {
      outline.style.opacity = "0";
      setTimeout(() => outline.remove(), 200);
    }, 500);
  }

  // Capture image and prepare for analysis
  async captureAndAnalyze(evt) {
    // Validate video element state
    if (!this.video) {
      throw new Error('Video element not found');
    }
    
    if (!this.currentStream) {
      throw new Error('No camera stream active');
    }
    
    if (this.video.videoWidth === 0 || this.video.videoHeight === 0) {
      throw new Error('Video stream not ready - please wait a moment and try again');
    }
    
    if (this.video.readyState < 2) {
      throw new Error('Video stream not ready - please wait a moment and try again');
    }
    
    const canvas = document.createElement("canvas");
    canvas.width = this.video.videoWidth;
    canvas.height = this.video.videoHeight;
    
    const ctx = canvas.getContext("2d");
    if (!ctx) {
      throw new Error('Failed to get canvas context');
    }
    
    try {
      ctx.drawImage(this.video, 0, 0, canvas.width, canvas.height);
    } catch (drawError) {
      throw new Error('Failed to capture image from camera - please try again');
    }
    
    const imageBase64 = canvas.toDataURL("image/jpeg", 0.9).split(",")[1];

    const scaleX = this.video.videoWidth / this.video.clientWidth;
    const scaleY = this.video.videoHeight / this.video.clientHeight;
    const clickX = Math.round(evt.clientX - this.video.getBoundingClientRect().left);
    const clickY = Math.round(evt.clientY - this.video.getBoundingClientRect().top);
    const actualX = Math.round(clickX * scaleX);
    const actualY = Math.round(clickY * scaleY);

    const crop = document.createElement("canvas");
    crop.width = crop.height = 100;
    const cropCtx = crop.getContext("2d");
    cropCtx.imageSmoothingEnabled = false;
    cropCtx.drawImage(canvas, actualX - 50, actualY - 50, 100, 100, 0, 0, 100, 100);

    const middleX = Math.floor(100 / 2);
    const middleY = Math.floor(100 / 2);
    const boxSize = Math.max(8, Math.floor(100 * 0.1));
    
    cropCtx.save();
    cropCtx.strokeStyle = "#ff0000";
    cropCtx.lineWidth = Math.max(2, Math.floor(100 * 0.02));
    cropCtx.strokeRect(middleX - boxSize / 2, middleY - boxSize / 2, boxSize, boxSize);
    cropCtx.restore();

    const croppedBase64 = crop.toDataURL("image/jpeg", 0.9).split(",")[1];

    return {
      imageBase64,
      croppedBase64,
      coordinates: { x: actualX, y: actualY, cropMiddleX: middleX, cropMiddleY: middleY, cropSize: 100 }
    };
  }

  // Setup camera switching functionality
  async setupSwitchCamera() {
    const devices = await navigator.mediaDevices.enumerateDevices();
    const videoDevices = devices.filter(d => d.kind === "videoinput");
    
    if (videoDevices.length > 1) {
      const switchCameraContainer = document.querySelector(".switch-camera-container");
      const switchCameraBtn = document.getElementById("switch-camera-btn");
      
      if (switchCameraContainer && switchCameraBtn) {
        switchCameraContainer.style.display = "block";
        switchCameraBtn.onclick = () => {
          this.facingMode = this.facingMode === "user" ? "environment" : "user";
          this.startCamera();
        };
      }
    }
  }

  // Update mobile reticle visibility
  updateReticleVisibility() {
    if (this.isMobile) {
      this.addMobileTargetingReticle();
    } else {
      this.removeMobileTargetingReticle();
    }
  }

  // Add mobile targeting reticle
  addMobileTargetingReticle() {
    this.removeMobileTargetingReticle();

    const reticle = document.createElement("div");
    reticle.className = "mobile-reticle";
    const template = document.getElementById("mobile-reticle-template");
    const clone = template.content.cloneNode(true);
    reticle.appendChild(clone);

    this.video.appendChild(reticle);

    setTimeout(() => {
      reticle.style.animation = "reticlePulse 2s ease-in-out infinite";
    }, 1000);
  }

  // Remove mobile targeting reticle
  removeMobileTargetingReticle() {
    const existingReticle = document.querySelector(".mobile-reticle");
    if (existingReticle) {
      existingReticle.remove();
    }
  }

  // Update camera instruction text
  updateCameraInstructions() {
    const instructionText = document.querySelector('.instruction-text');
    if (!instructionText) return;
    
    if (this.isCameraReady()) {
      instructionText.textContent = 'Center object in circle & tap, or type name above';
      instructionText.style.color = '#ffffff';
    } else if (this.currentStream) {
      instructionText.textContent = 'Camera initializing... please wait';
      instructionText.style.color = '#ffa500';
    } else if (this.isRequestingPermission) {
      instructionText.textContent = 'Requesting camera permission...';
      instructionText.style.color = '#ffa500';
    } else {
      instructionText.textContent = 'Camera not available - check permissions';
      instructionText.style.color = '#ff6b6b';
    }
  }

  // Handle camera errors
  handleCameraError(err) {
    let errorMessage;
    
    if (err.name === "NotAllowedError") {
      errorMessage = "Camera access denied. Please allow camera access in browser settings.";
    } else if (err.name === "NotFoundError") {
      errorMessage = "No camera found on this device.";
    } else if (this.isSafari && err.name === "NotSupportedError") {
      errorMessage = "Camera not supported in Safari. Try using Chrome or Firefox.";
    } else {
      errorMessage = `Camera error: ${err.message}`;
    }
    
    this.showError(errorMessage);
  }

  // Show error message
  showError(message) {
    const msgBox = document.querySelector(".permission-message");
    if (msgBox) {
      msgBox.hidden = false;
      msgBox.textContent = message;
    }
  }

  // Hide error message
  hideError() {
    const msgBox = document.querySelector(".permission-message");
    if (msgBox) {
      msgBox.hidden = true;
    }
  }

  // Create closable error message
  createClosableErrorMessage(message) {
    const snapshots = document.querySelector(".snapshots-container");
    const errorDiv = uiManager.createErrorMessage(message, snapshots);
    return errorDiv;
  }

  // Request camera permission with persistence
  async requestPermission() {
    this.isRequestingPermission = true;
    this.updateCameraInstructions();
    
    try {
      const stream = await navigator.mediaDevices.getUserMedia({ video: true });
      this.saveCameraPermission(true);
      stream.getTracks().forEach(track => track.stop());
      this.isRequestingPermission = false;
      await this.startCamera();
      return true;
    } catch (error) {
      this.isRequestingPermission = false;
      this.saveCameraPermission(false);
      this.updateCameraInstructions();
      return false;
    }
  }

  // Check payment setup for sidebar-based system
  async checkPaymentSetupForSidebar() {
    // Check if dev mode is enabled (localhost auto-enable)
    if (window.app && window.app.hasPaymentSetup === true) {
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

  // Increment usage for sidebar-based system
  async incrementUsageForSidebar() {
    const deviceToken = localStorage.getItem('molDeviceToken');
    if (!deviceToken) return;
    
    try {
      const response = await fetch('/increment-usage', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ device_token: deviceToken })
      });
      
      if (response.ok) {
        const result = await response.json();
        const cardInfo = JSON.parse(localStorage.getItem('molCardInfo') || '{}');
        cardInfo.usage = result.usage;
        localStorage.setItem('molCardInfo', JSON.stringify(cardInfo));
      }
    } catch (error) {
      logger.warn('Usage increment failed', error);
    }
  }

  // Cleanup resources
  cleanup() {
    this.stopCurrentStream();
    this.removeMobileTargetingReticle();
  }
}

// Export singleton instance
export const cameraManager = new CameraManager(); 