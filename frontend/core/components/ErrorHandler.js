// ErrorHandler.js - Centralized error handling and display
export class ErrorHandler {
  constructor() {
    this.errorPanel = null;
    this.errorList = null;
    this.errors = [];
  }

  initialize() {
    this.errorPanel = document.getElementById('error-panel');
    this.errorList = document.getElementById('error-list');
    
    if (this.errorPanel && this.errorList) {
      this.setupEventListeners();
      console.log('✅ Error handler initialized');
    } else {
      console.warn('⚠️ Error panel elements not found');
    }
  }

  setupEventListeners() {
    // Listen for app errors
    document.addEventListener('appError', (e) => {
      this.handleError(e.detail.message, e.detail.type);
    });

    // Clear errors button
    const clearBtn = document.getElementById('clear-errors');
    if (clearBtn) {
      clearBtn.addEventListener('click', () => this.clearErrors());
    }

    // Toggle errors button
    const toggleBtn = document.getElementById('toggle-errors');
    if (toggleBtn) {
      toggleBtn.addEventListener('click', () => this.toggleErrorPanel());
    }
  }

  handleError(message, type = 'error') {
    const error = {
      id: Date.now(),
      message,
      type,
      timestamp: new Date().toISOString()
    };

    this.errors.push(error);
    this.displayError(error);
    this.showErrorPanel();

    // Log to console
    console.error(`❌ [${type.toUpperCase()}] ${message}`);

    // Auto-hide after 10 seconds for non-critical errors
    if (type !== 'critical') {
      setTimeout(() => {
        this.removeError(error.id);
      }, 10000);
    }
  }

  displayError(error) {
    if (!this.errorList) return;

    const errorItem = document.createElement('div');
    errorItem.className = `error-item ${error.type}`;
    errorItem.dataset.errorId = error.id;
    errorItem.style.cssText = `
      padding: 10px;
      margin-bottom: 8px;
      border-radius: 4px;
      font-size: 14px;
      line-height: 1.4;
      display: flex;
      justify-content: space-between;
      align-items: flex-start;
      gap: 10px;
    `;

    // Set background color based on error type
    switch (error.type) {
      case 'critical':
        errorItem.style.background = 'rgba(255, 68, 68, 0.2)';
        errorItem.style.border = '1px solid #ff4444';
        break;
      case 'warning':
        errorItem.style.background = 'rgba(255, 170, 0, 0.2)';
        errorItem.style.border = '1px solid #ffaa00';
        break;
      default:
        errorItem.style.background = 'rgba(255, 68, 68, 0.1)';
        errorItem.style.border = '1px solid #ff6666';
    }

    const messageDiv = document.createElement('div');
    messageDiv.textContent = error.message;
    messageDiv.style.flex = '1';

    const closeBtn = document.createElement('button');
    closeBtn.innerHTML = '✕';
    closeBtn.style.cssText = `
      background: none;
      border: none;
      color: inherit;
      cursor: pointer;
      font-size: 16px;
      padding: 2px;
      opacity: 0.7;
      transition: opacity 0.2s;
    `;
    closeBtn.onmouseover = () => closeBtn.style.opacity = '1';
    closeBtn.onmouseout = () => closeBtn.style.opacity = '0.7';
    closeBtn.onclick = () => this.removeError(error.id);

    errorItem.appendChild(messageDiv);
    errorItem.appendChild(closeBtn);
    this.errorList.appendChild(errorItem);
  }

  removeError(errorId) {
    // Remove from errors array
    this.errors = this.errors.filter(e => e.id !== errorId);

    // Remove from DOM
    const errorElement = this.errorList?.querySelector(`[data-error-id="${errorId}"]`);
    if (errorElement) {
      errorElement.remove();
    }

    // Hide panel if no errors left
    if (this.errors.length === 0) {
      this.hideErrorPanel();
    }
  }

  clearErrors() {
    this.errors = [];
    if (this.errorList) {
      this.errorList.innerHTML = '';
    }
    this.hideErrorPanel();
  }

  showErrorPanel() {
    if (this.errorPanel) {
      this.errorPanel.classList.remove('hidden');
    }
  }

  hideErrorPanel() {
    if (this.errorPanel) {
      this.errorPanel.classList.add('hidden');
    }
  }

  toggleErrorPanel() {
    if (this.errorPanel) {
      this.errorPanel.classList.toggle('hidden');
      
      const toggleBtn = document.getElementById('toggle-errors');
      if (toggleBtn) {
        toggleBtn.textContent = this.errorPanel.classList.contains('hidden') ? 'Show' : 'Hide';
      }
    }
  }

  // Global error handler for uncaught errors
  setupGlobalErrorHandler() {
    window.addEventListener('error', (event) => {
      this.handleError(`Uncaught error: ${event.message}`, 'critical');
    });

    window.addEventListener('unhandledrejection', (event) => {
      this.handleError(`Unhandled promise rejection: ${event.reason}`, 'critical');
    });
  }

  // Get error statistics
  getErrorStats() {
    const stats = {
      total: this.errors.length,
      critical: this.errors.filter(e => e.type === 'critical').length,
      warning: this.errors.filter(e => e.type === 'warning').length,
      error: this.errors.filter(e => e.type === 'error').length
    };
    return stats;
  }
} 