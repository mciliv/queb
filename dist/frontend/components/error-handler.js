// Error Handler - Centralized AI & Human Readable Error Reporting System
// Auto-generated based on error-reporting.mdc specifications

class ErrorReporter {
  constructor() {
    this.errors = [];
    this.maxErrors = 50;
    this.init();
  }

  init() {
    // Capture all JavaScript errors
    window.addEventListener('error', (event) => {
      this.reportError({
        type: 'JavaScript Error',
        message: event.message,
        source: `${event.filename}:${event.lineno}:${event.colno}`,
        stack: event.error?.stack,
        timestamp: new Date()
      });
    });

    // Capture unhandled promise rejections
    window.addEventListener('unhandledrejection', (event) => {
      this.reportError({
        type: 'Unhandled Promise Rejection',
        message: event.reason?.message || event.reason,
        source: 'Promise',
        stack: event.reason?.stack,
        timestamp: new Date()
      });
    });

    // Setup UI handlers
    this.setupUI();
  }

  reportError(error) {
    // Add to errors array
    this.errors.unshift(error);
    
    // Limit array size
    if (this.errors.length > this.maxErrors) {
      this.errors = this.errors.slice(0, this.maxErrors);
    }

    // Update UI
    this.updateErrorPanel();

    // Log to console for AI parsing - STRUCTURED JSON FORMAT
    console.error('🚨 ERROR REPORTED:', {
      timestamp: error.timestamp.toISOString(),
      type: error.type,
      message: error.message,
      source: error.source,
      stack: error.stack
    });

    // Send to server for logging (optional)
    this.sendToServer(error);
  }

  // Manual error reporting for custom errors
  report(type, message, source = 'Manual', stack = null) {
    this.reportError({
      type,
      message,
      source,
      stack,
      timestamp: new Date()
    });
  }

  updateErrorPanel() {
    const panel = document.getElementById('error-panel');
    const errorList = document.getElementById('error-list');
    
    if (!panel || !errorList) return;

    // Show panel if hidden and has errors
    if (this.errors.length > 0) {
      panel.classList.remove('hidden');
    }

    // Update error list
    errorList.innerHTML = this.errors.map(error => `
      <div class="error-item">
        <div class="error-timestamp">${error.timestamp.toLocaleString()}</div>
        <div class="error-type">${error.type}</div>
        <div class="error-message">${this.escapeHtml(error.message)}</div>
        <div class="error-source">${this.escapeHtml(error.source)}</div>
        ${error.stack ? `<div class="error-stack">${this.escapeHtml(error.stack)}</div>` : ''}
      </div>
    `).join('');
  }

  escapeHtml(text) {
    const div = document.createElement('div');
    div.textContent = text;
    return div.innerHTML;
  }

  setupUI() {
    // Wait for DOM to be ready
    if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', () => this.setupUIHandlers());
    } else {
      this.setupUIHandlers();
    }
  }

  setupUIHandlers() {
    // Clear errors button
    document.addEventListener('click', (e) => {
      if (e.target.id === 'clear-errors') {
        this.clearErrors();
      }
      if (e.target.id === 'toggle-errors') {
        this.togglePanel();
      }
    });
  }

  clearErrors() {
    this.errors = [];
    this.updateErrorPanel();
    const panel = document.getElementById('error-panel');
    if (panel) {
      panel.classList.add('hidden');
    }
  }

  togglePanel() {
    const panel = document.getElementById('error-panel');
    if (panel) {
      panel.classList.toggle('hidden');
    }
  }

  sendToServer(error) {
    // Optional: Send errors to backend for logging
    fetch('/api/log-error', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(error)
    }).catch(() => {
      // Fail silently if server unavailable
    });
  }

  // AI-readable error summary
  getErrorSummary() {
    return {
      totalErrors: this.errors.length,
      recentErrors: this.errors.slice(0, 5),
      errorTypes: [...new Set(this.errors.map(e => e.type))],
      commonSources: this.getMostCommonSources()
    };
  }

  getMostCommonSources() {
    const sources = this.errors.map(e => e.source);
    const counts = {};
    sources.forEach(source => counts[source] = (counts[source] || 0) + 1);
    return Object.entries(counts)
      .sort(([,a], [,b]) => b - a)
      .slice(0, 5)
      .map(([source, count]) => ({ source, count }));
  }

  // Test error reporting (for debugging)
  testError() {
    this.report('Test Error', 'This is a test error for debugging', 'error-handler.js');
  }
}

// Initialize global error reporter
window.errorReporter = new ErrorReporter();

// Convenience function for manual error reporting
window.reportError = (type, message, source) => {
  window.errorReporter.report(type, message, source);
};

// Add to global debug functions
window.testError = () => window.errorReporter.testError();


 to report custom errors');
 for AI analysis');
 to test the error panel'); 