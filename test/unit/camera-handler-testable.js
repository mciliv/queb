// Test-specific version of camera handler with injected dependencies
class CameraHandler {
  constructor(uiManager, paymentManager) {
    this.uiManager = uiManager;
    this.paymentManager = paymentManager;
    this.isMobile = /Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini/i.test(navigator.userAgent);
  }

  // Handle photo upload analysis
  async handlePhotoUpload(e) {
    const file = e.target.files[0];
    if (!file) return;

    if (!file.type.startsWith("image/")) {
      alert("Please select an image file");
      return;
    }

    await this.displayUploadedImage(file);
  }

  // Handle URL-based image analysis
  async handleUrlAnalysis() {
    const photoUrl = document.getElementById("photo-url");
    const url = photoUrl.value.trim();
    
    if (!url) {
      alert("Please enter an image URL");
      return;
    }

    try {
      new URL(url);
    } catch {
      alert("Please enter a valid URL");
      return;
    }

    try {
      const imageBase64 = await this.uiManager.urlToBase64(url);
      const byteCharacters = atob(imageBase64);
      const byteNumbers = new Array(byteCharacters.length);
      
      for (let i = 0; i < byteCharacters.length; i++) {
        byteNumbers[i] = byteCharacters.charCodeAt(i);
      }
      
      const byteArray = new Uint8Array(byteNumbers);
      const blob = new Blob([byteArray], { type: "image/jpeg" });
      const file = new File([blob], "url-image.jpg", { type: "image/jpeg" });

      await this.displayUploadedImage(file);
      photoUrl.value = "";
      
    } catch (err) {
      this.createClosableErrorMessage(`Error loading image from URL: ${err.message}`);
    }
  }

  // Display uploaded image for interactive analysis
  async displayUploadedImage(file) {
    const photoOptions = document.getElementById("photo-options");
    photoOptions.innerHTML = "";

    const imageContainer = document.createElement("div");
    imageContainer.className = "uploaded-image-container";

    const img = document.createElement("img");
    
    try {
      const imageBase64 = await this.uiManager.fileToBase64(file);
      img.src = `data:${file.type};base64,${imageBase64}`;
      img.dataset.imageBase64 = imageBase64;
      
      if (this.isMobile) {
        const crosshair = document.createElement("div");
        crosshair.className = "mobile-crosshair";
        imageContainer.appendChild(crosshair);
      }
      
      imageContainer.appendChild(img);
      photoOptions.appendChild(imageContainer);
      
      img.addEventListener("click", (e) => this.handleImageClick(e, img));
    } catch (error) {
      console.error("Error processing image:", error.message);
    }
  }

  // Handle image click for molecular analysis
  async handleImageClick(evt, img) {
    try {
      // Check payment method first
      try {
        await this.paymentManager.checkPaymentMethod();
      } catch (paymentError) {
        console.warn("Payment check failed, proceeding anyway:", paymentError);
      }

      if (!img.dataset.imageBase64) {
        this.createClosableErrorMessage("No image data available");
        return;
      }
      
      const rect = img.getBoundingClientRect();
      const clickX = evt.clientX - rect.left;
      const clickY = evt.clientY - rect.top;
      const relativeX = clickX / rect.width;
      const relativeY = clickY / rect.height;

      // Create temporary image to get natural dimensions
      const tempImg = new Image();
      tempImg.src = img.src;
      
      await new Promise((resolve, reject) => {
        tempImg.onload = resolve;
        tempImg.onerror = () => reject(new Error("Failed to process image"));
      });

      const x = Math.round(relativeX * tempImg.width);
      const y = Math.round(relativeY * tempImg.height);

      const response = await fetch('/image-molecules', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          imageBase64: img.dataset.imageBase64,
          x: x,
          y: y
        })
      });

      if (!response.ok) {
        throw new Error(`HTTP ${response.status}: ${response.statusText}`);
      }

      const result = await response.json();
      
      this.emitAnalysisResult(result, 'Photo', 'Unknown', false, img.dataset.imageBase64);
      
      try {
        await this.paymentManager.incrementUsage();
      } catch (error) {
        console.warn("Usage increment failed:", error);
      }
      
    } catch (error) {
      this.createClosableErrorMessage(`Error: ${error.message}`);
    }
  }

  // Emit analysis result as custom event
  emitAnalysisResult(output, icon = 'Photo', objectName = 'Unknown', useQuotes = false, croppedImageData = null) {
    const event = new CustomEvent('imageAnalysisComplete', {
      detail: {
        output,
        icon,
        objectName,
        useQuotes,
        croppedImageData
      }
    });
    
    document.dispatchEvent(event);
  }

  // Create closable error message
  createClosableErrorMessage(message) {
    const snapshotsContainer = document.querySelector('.snapshots-container');
    return this.uiManager.createErrorMessage(message, snapshotsContainer);
  }

  // Update scroll handles
  updateScrollHandles() {
    if (typeof window.updateScrollHandles === 'function') {
      window.updateScrollHandles();
    }
  }

  // Setup event listeners
  setupEventListeners() {
    const photoUpload = document.getElementById("photo-upload");
    if (photoUpload) {
      photoUpload.addEventListener("change", (e) => this.handlePhotoUpload(e));
    }

    const photoUrl = document.getElementById("photo-url");
    const urlAnalyze = document.getElementById("url-analyze");
    if (photoUrl && urlAnalyze) {
      photoUrl.addEventListener("keyup", (e) => {
        if (e.key === "Enter") this.handleUrlAnalysis();
      });
      urlAnalyze.addEventListener("click", () => this.handleUrlAnalysis());
    }
  }
}

module.exports = { CameraHandler };