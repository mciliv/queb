/**
 * Error Screenshot Utility - Frontend screenshot capture for debugging
 *
 * This utility automatically captures screenshots when JavaScript errors occur,
 * helping with debugging and issue reproduction in the Cursor development environment.
 *
 * INTEGRATION:
 * - Automatically hooks into window.onerror
 * - Captures screenshots using html2canvas or browser APIs
 * - Sends screenshots to backend via existing error logging endpoint
 * - Only active when ENABLE_ERROR_SCREENSHOTS is enabled
 */

class ErrorScreenshot {
  constructor(options = {}) {
    this.apiEndpoint = options.apiEndpoint || '/api/log-error';
    this.enabled = options.enabled || false;
    this.maxRetries = options.maxRetries || 2;
    this.retryDelay = options.retryDelay || 1000;

    // Initialize if enabled
    if (this.enabled) {
      this.init();
    }
  }

  /**
   * Initialize screenshot capturing
   */
  init() {
    if (!this.enabled) return;

    // Hook into global error handler
    this.originalOnError = window.onerror;
    window.onerror = this.handleError.bind(this);

    // Hook into unhandled promise rejections
    window.addEventListener('unhandledrejection', this.handleUnhandledRejection.bind(this));

    console.log('🔧 Error screenshot capture enabled');
  }

  /**
   * Handle JavaScript errors
   */
  async handleError(message, source, lineno, colno, error) {
    try {
      // Call original error handler if it exists
      if (this.originalOnError) {
        this.originalOnError(message, source, lineno, colno, error);
      }

      // Capture screenshot and send error
      await this.captureAndSendError({
        message: message,
        source: source,
        line: lineno,
        column: colno,
        stack: error?.stack,
        type: 'javascript_error',
        url: window.location.href,
        userAgent: navigator.userAgent,
        timestamp: new Date().toISOString()
      });

    } catch (screenshotError) {
      console.error('Failed to capture error screenshot:', screenshotError);
    }
  }

  /**
   * Handle unhandled promise rejections
   */
  async handleUnhandledRejection(event) {
    try {
      const error = event.reason;
      await this.captureAndSendError({
        message: error?.message || 'Unhandled Promise Rejection',
        stack: error?.stack,
        type: 'unhandled_rejection',
        reason: error,
        url: window.location.href,
        userAgent: navigator.userAgent,
        timestamp: new Date().toISOString()
      });

    } catch (screenshotError) {
      console.error('Failed to capture rejection screenshot:', screenshotError);
    }
  }

  /**
   * Capture screenshot and send error data
   */
  async captureAndSendError(errorInfo) {
    if (!this.enabled) return;

    try {
      // Capture screenshot
      const screenshotData = await this.captureScreenshot();

      if (screenshotData) {
        // Send error with screenshot to backend
        await this.sendErrorWithScreenshot(errorInfo, screenshotData);
      } else {
        // Send error without screenshot if capture failed
        await this.sendError(errorInfo);
      }

    } catch (error) {
      console.error('Error screenshot capture failed:', error);
      // Still try to send basic error info
      try {
        await this.sendError(errorInfo);
      } catch (sendError) {
        console.error('Failed to send error data:', sendError);
      }
    }
  }

  /**
   * Capture screenshot of current page
   */
  async captureScreenshot() {
    // Try html2canvas first (if available)
    if (typeof html2canvas !== 'undefined') {
      try {
        const canvas = await html2canvas(document.body, {
          width: window.innerWidth,
          height: window.innerHeight,
          useCORS: true,
          allowTaint: false,
          scale: 0.8 // Reduce quality for performance
        });

        return canvas.toDataURL('image/png');
      } catch (error) {
        console.warn('html2canvas capture failed:', error);
      }
    }

    // Fallback: Try browser native screenshot API (limited browser support)
    if (navigator.mediaDevices && navigator.mediaDevices.getDisplayMedia) {
      try {
        console.log('Attempting native screenshot API...');
        // Note: This would require user permission and is complex to implement
        // For now, we'll skip this fallback
      } catch (error) {
        console.warn('Native screenshot API failed:', error);
      }
    }

    console.warn('No screenshot method available');
    return null;
  }

  /**
   * Send error data with screenshot to backend
   */
  async sendErrorWithScreenshot(errorInfo, screenshotData) {
    const payload = {
      ...errorInfo,
      screenshot: {
        data: screenshotData,
        format: 'png',
        capturedAt: new Date().toISOString()
      }
    };

    await this.sendError(payload);
  }

  /**
   * Send error data to backend with retry logic
   */
  async sendError(errorData, attempt = 1) {
    try {
      const response = await fetch(this.apiEndpoint, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(errorData)
      });

      if (!response.ok) {
        throw new Error(`HTTP ${response.status}: ${response.statusText}`);
      }

      const result = await response.json();
      console.log('Error data sent successfully:', result);

    } catch (error) {
      console.error(`Error send attempt ${attempt} failed:`, error);

      // Retry logic
      if (attempt < this.maxRetries) {
        console.log(`Retrying in ${this.retryDelay}ms...`);
        await new Promise(resolve => setTimeout(resolve, this.retryDelay));
        return this.sendError(errorData, attempt + 1);
      }

      throw error;
    }
  }

  /**
   * Enable screenshot capturing
   */
  enable() {
    this.enabled = true;
    this.init();
  }

  /**
   * Disable screenshot capturing
   */
  disable() {
    this.enabled = false;

    // Restore original error handler
    if (this.originalOnError) {
      window.onerror = this.originalOnError;
    }

    // Remove unhandled rejection listener
    window.removeEventListener('unhandledrejection', this.handleUnhandledRejection);
  }

  /**
   * Manually trigger screenshot capture (for testing)
   */
  async testScreenshot() {
    try {
      console.log('Testing screenshot capture...');
      const screenshot = await this.captureScreenshot();

      if (screenshot) {
        console.log('Screenshot captured successfully');
        console.log('Screenshot preview:', screenshot.substring(0, 50) + '...');

        // Test sending
        await this.captureAndSendError({
          message: 'Test screenshot capture',
          type: 'test_screenshot',
          url: window.location.href,
          timestamp: new Date().toISOString()
        });

        return true;
      } else {
        console.warn('Screenshot capture failed');
        return false;
      }

    } catch (error) {
      console.error('Screenshot test failed:', error);
      return false;
    }
  }
}

// Export for use in other scripts
window.ErrorScreenshot = ErrorScreenshot;