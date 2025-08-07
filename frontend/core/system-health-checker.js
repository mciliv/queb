class AppChecker {
  constructor() {
    this.results = [];
    this.errors = [];
  }

  logResult(test, status, message, data = null) {
    const result = {
      test,
      status,
      message,
      timestamp: new Date().toISOString(),
      data
    };
    
    this.results.push(result);
    console.log(`[${status}] ${test}: ${message}`);
    
    if (data) {
      console.log('Data:', data);
    }
    
    return result;
  }

  async checkAllComponents() {
    console.log('🔍 Starting comprehensive app check...');
    
    this.checkDOMElements();
    this.checkModules();
    this.checkPaymentSystem();
    this.checkCameraSystem();
    this.checkEventListeners();
    this.checkAPIEndpoints();
    this.checkUIComponents();
    this.checkStorage();
    this.checkPermissions();
    
    this.generateReport();
  }

  checkDOMElements() {
    const elements = [
      'video-feed',
      'object-input',
      'snapshots-container',
      'payment-modal',

      'account-status',
      'account-name'
    ];
    
    const missing = [];
    const found = [];
    
    elements.forEach(id => {
      const element = document.getElementById(id);
      if (element) {
        found.push(id);
      } else {
        missing.push(id);
      }
    });
    
    if (missing.length === 0) {
      this.logResult('DOM Elements', 'PASS', 'All required elements found', { found });
    } else {
      this.logResult('DOM Elements', 'FAIL', `Missing elements: ${missing.join(', ')}`, { found, missing });
    }
  }

  checkModules() {
    const modules = [
      'paymentManager',
      'cameraManager', 
      'cameraHandler',
      'uiManager'
    ];
    
    const missing = [];
    const found = [];
    
    modules.forEach(module => {
      if (window[module]) {
        found.push(module);
      } else {
        missing.push(module);
      }
    });
    
    if (missing.length === 0) {
      this.logResult('Modules', 'PASS', 'All modules loaded', { found });
    } else {
      this.logResult('Modules', 'FAIL', `Missing modules: ${missing.join(', ')}`, { found, missing });
    }
  }

  checkPaymentSystem() {
    if (!window.paymentManager) {
      this.logResult('Payment System', 'FAIL', 'Payment manager not available');
      return;
    }
    
    const checks = [
      { name: 'isDeveloperAccount', test: () => typeof window.paymentManager.isDeveloperAccount === 'function' },
      { name: 'checkPaymentMethod', test: () => typeof window.paymentManager.checkPaymentMethod === 'function' },
      { name: 'setupDeveloperAccount', test: () => typeof window.paymentManager.setupDeveloperAccount === 'function' }
    ];
    
    const passed = [];
    const failed = [];
    
    checks.forEach(check => {
      if (check.test()) {
        passed.push(check.name);
      } else {
        failed.push(check.name);
      }
    });
    
    if (failed.length === 0) {
      this.logResult('Payment System', 'PASS', 'All payment functions available', { passed });
    } else {
      this.logResult('Payment System', 'FAIL', `Missing functions: ${failed.join(', ')}`, { passed, failed });
    }
  }

  checkCameraSystem() {
    if (!window.cameraManager) {
      this.logResult('Camera System', 'FAIL', 'Camera manager not available');
      return;
    }
    
    const checks = [
      { name: 'initialize', test: () => typeof window.cameraManager.initialize === 'function' },
      { name: 'startCamera', test: () => typeof window.cameraManager.startCamera === 'function' },
      { name: 'isCameraReady', test: () => typeof window.cameraManager.isCameraReady === 'function' }
    ];
    
    const passed = [];
    const failed = [];
    
    checks.forEach(check => {
      if (check.test()) {
        passed.push(check.name);
      } else {
        failed.push(check.name);
      }
    });
    
    if (failed.length === 0) {
      this.logResult('Camera System', 'PASS', 'All camera functions available', { passed });
    } else {
      this.logResult('Camera System', 'FAIL', `Missing functions: ${failed.join(', ')}`, { passed, failed });
    }
  }

  checkEventListeners() {
    const video = document.getElementById('video-feed');
    const input = document.getElementById('object-input');
    
    const checks = [];
    
    if (video) {
      checks.push({ name: 'Video click listener', element: video, event: 'click' });
      checks.push({ name: 'Video touch listener', element: video, event: 'touchstart' });
    }
    
    if (input) {
      checks.push({ name: 'Input keyup listener', element: input, event: 'keyup' });
      checks.push({ name: 'Input focus listener', element: input, event: 'focus' });
    }
    
    const passed = [];
    const failed = [];
    
    checks.forEach(check => {
      const listeners = getEventListeners(check.element, check.event);
      if (listeners && listeners.length > 0) {
        passed.push(check.name);
      } else {
        failed.push(check.name);
      }
    });
    
    if (failed.length === 0) {
      this.logResult('Event Listeners', 'PASS', 'All event listeners attached', { passed });
    } else {
      this.logResult('Event Listeners', 'WARN', `Missing listeners: ${failed.join(', ')}`, { passed, failed });
    }
  }

  async checkAPIEndpoints() {
    const endpoints = [
      '/text-molecules',
      '/image-molecules',
      '/validate-payment',
      '/increment-usage'
    ];
    
    const results = [];
    
    for (const endpoint of endpoints) {
      try {
        const response = await fetch(endpoint, { method: 'HEAD' });
        results.push({ endpoint, status: response.status, available: response.status !== 404 });
      } catch (error) {
        results.push({ endpoint, status: 'error', available: false, error: error.message });
      }
    }
    
    const available = results.filter(r => r.available);
    const unavailable = results.filter(r => !r.available);
    
    if (unavailable.length === 0) {
      this.logResult('API Endpoints', 'PASS', 'All endpoints available', { available });
    } else {
      this.logResult('API Endpoints', 'FAIL', `Unavailable endpoints: ${unavailable.map(r => r.endpoint).join(', ')}`, { available, unavailable });
    }
  }

  checkUIComponents() {
    const components = [
      { name: 'Payment Popdown', selector: '#payment-modal', shouldExist: true },

      { name: 'Camera Container', selector: '#camera-container', shouldExist: true },
      { name: 'Photo Options', selector: '#photo-options', shouldExist: true }
    ];
    
    const passed = [];
    const failed = [];
    
    components.forEach(component => {
      const element = document.querySelector(component.selector);
      if (element && component.shouldExist) {
        passed.push(component.name);
      } else if (!element && !component.shouldExist) {
        passed.push(component.name);
      } else {
        failed.push(component.name);
      }
    });
    
    if (failed.length === 0) {
      this.logResult('UI Components', 'PASS', 'All UI components present', { passed });
    } else {
      this.logResult('UI Components', 'FAIL', `Missing components: ${failed.join(', ')}`, { passed, failed });
    }
  }

  checkStorage() {
    const storageKeys = [
      'molDeviceToken',
      'molCardInfo',
      'molDeveloperUser',
      'mol_device_id',
      'mol_camera_permission'
    ];
    
    const found = [];
    const missing = [];
    
    storageKeys.forEach(key => {
      const value = localStorage.getItem(key);
      if (value) {
        found.push(key);
      } else {
        missing.push(key);
      }
    });
    
    if (missing.length === 0) {
      this.logResult('Local Storage', 'PASS', 'All storage keys present', { found });
    } else {
      this.logResult('Local Storage', 'WARN', `Missing keys: ${missing.join(', ')}`, { found, missing });
    }
  }

  checkPermissions() {
    const permissions = [
      { name: 'Camera', api: 'camera' },
      { name: 'Media Devices', api: 'mediaDevices' }
    ];
    
    const available = [];
    const unavailable = [];
    
    permissions.forEach(permission => {
      if (navigator.permissions && navigator.permissions.query) {
        available.push(permission.name);
      } else {
        unavailable.push(permission.name);
      }
    });
    
    if (unavailable.length === 0) {
      this.logResult('Permissions', 'PASS', 'All permissions available', { available });
    } else {
      this.logResult('Permissions', 'WARN', `Unavailable permissions: ${unavailable.join(', ')}`, { available, unavailable });
    }
  }

  generateReport() {
    const total = this.results.length;
    const passed = this.results.filter(r => r.status === 'PASS').length;
    const failed = this.results.filter(r => r.status === 'FAIL').length;
    const warnings = this.results.filter(r => r.status === 'WARN').length;
    
    console.log('\n📊 App Check Report');
    console.log('==================');
    console.log(`Total Tests: ${total}`);
    console.log(`✅ Passed: ${passed}`);
    console.log(`❌ Failed: ${failed}`);
    console.log(`⚠️ Warnings: ${warnings}`);
    console.log(`Success Rate: ${((passed / total) * 100).toFixed(1)}%`);
    
    if (failed > 0) {
      console.log('\n❌ Failed Tests:');
      this.results.filter(r => r.status === 'FAIL').forEach(result => {
        console.log(`  - ${result.test}: ${result.message}`);
      });
    }
    
    if (warnings > 0) {
      console.log('\n⚠️ Warnings:');
      this.results.filter(r => r.status === 'WARN').forEach(result => {
        console.log(`  - ${result.test}: ${result.message}`);
      });
    }
    
    console.log('\n📋 Full Results:', this.results);
  }
}

window.appChecker = new AppChecker();

console.log('✅ App Checker loaded. Run: appChecker.checkAllComponents()'); 