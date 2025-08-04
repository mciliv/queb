// AppShell.js - Main application shell and keyboard shortcuts
import { InputSection } from './InputSection.js';

export class AppShell {
  constructor() {
    this.inputSection = new InputSection();
    this.snapshots = null;
    this.isProcessing = false;
    this.hasPaymentSetup = false;
    this.paymentEnabled = false;
    this.currentAnalysisType = null;
    this.lastAnalysis = null;
  }

  async initialize() {
    this.snapshots = document.querySelector(".snapshots-container");
    
    // Initialize input section
    this.inputSection.initialize();

    this.setupEventListeners();
    this.setupKeyboardShortcuts();
    this.setupKeyboardHint();

    console.log('✅ App shell initialized');
  }

  setupEventListeners() {
    // Text analysis - handled by InputSection
    this.inputSection.objectInput.addEventListener("keyup", async (e) => {
      if (e.key !== "Enter") return;
      await this.handleTextAnalysis();
    });

    // Image analysis events
    document.addEventListener('imageAnalysisComplete', (e) => {
      const { output, icon, objectName, useQuotes, croppedImageData } = e.detail;
      this.processAnalysisResult(output, icon, objectName, useQuotes, croppedImageData);
    });
  }

  setupKeyboardShortcuts() {
    // Keyboard shortcut (Cmd+K / Ctrl+K)
    document.addEventListener('keydown', (event) => {
      if ((event.metaKey || event.ctrlKey) && event.key === 'k') {
        event.preventDefault();
        this.inputSection.focusInput();
      }
    });
  }

  setupKeyboardHint() {
    const hintKey = document.getElementById('hint-key');
    if (!hintKey) return;

    const isMac = navigator.platform.toUpperCase().indexOf('MAC') >= 0;
    hintKey.textContent = isMac ? '⌘K' : 'Ctrl+K';

    // Show hint only on non-touch devices
    if ('ontouchstart' in window) {
      const keyboardHint = document.getElementById('keyboard-hint');
      if (keyboardHint) {
        keyboardHint.style.display = 'none';
      }
    }
  }

  async handleTextAnalysis() {
    if (this.isProcessing) return;

    const inputValue = this.inputSection.getInputValue();
    if (!inputValue) return;

    this.showProcessing();
    this.currentAnalysisType = 'text';

    try {
      const response = await fetch('/api/analyze-text', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ text: inputValue })
      });

      if (!response.ok) throw new Error('Analysis failed');
      
      const result = await response.json();
      this.processAnalysisResult(result.output, null, inputValue, false, null);
    } catch (error) {
      this.showError('Text analysis failed: ' + error.message);
    } finally {
      this.hideProcessing();
    }
  }

  showProcessing() {
    this.isProcessing = true;
    this.inputSection.disableInput();
  }

  hideProcessing() {
    this.isProcessing = false;
    this.inputSection.enableInput();
  }

  showError(message) {
    console.error('❌ Error:', message);
    // Dispatch error event for UI components to handle
    document.dispatchEvent(new CustomEvent('appError', { 
      detail: { message, type: 'error' } 
    }));
  }

  async processAnalysisResult(output, icon, objectName, useQuotes, croppedImageData) {
    // Dispatch event for other components to handle
    document.dispatchEvent(new CustomEvent('analysisResult', {
      detail: { output, icon, objectName, useQuotes, croppedImageData }
    }));
  }

  clearResults() {
    if (this.snapshots) {
      this.snapshots.innerHTML = '';
    }
    const gldiv = document.getElementById('gldiv');
    if (gldiv) {
      gldiv.innerHTML = '';
    }
  }

  // Payment management
  setPaymentStatus(enabled, hasSetup) {
    this.paymentEnabled = enabled;
    this.hasPaymentSetup = hasSetup;
  }

  hidePaymentSection() {
    const paymentSection = document.getElementById('payment-section');
    if (paymentSection) {
      paymentSection.style.display = 'none';
    }
  }

  showPaymentSection() {
    const paymentSection = document.getElementById('payment-section');
    if (paymentSection) {
      paymentSection.style.display = 'block';
    }
  }
} 