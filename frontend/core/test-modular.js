// test-modular.js - Test script for modular components
console.log('🧪 Testing modular components...');

// Wait for app to initialize
setTimeout(() => {
  if (window.molecularApp) {
    console.log('✅ Modular app initialized');
    
    // Test component access
    const appShell = window.molecularApp.getAppShell();
    const molecularViewer = window.molecularApp.getMolecularViewer();
    const errorHandler = window.molecularApp.getErrorHandler();
    
    console.log('✅ Components accessible:', {
      appShell: !!appShell,
      molecularViewer: !!molecularViewer,
      errorHandler: !!errorHandler
    });
    
    // Test keyboard shortcut simulation
    console.log('⌨️ Testing keyboard shortcut (Cmd+K)...');
    const cmdKEvent = new KeyboardEvent('keydown', {
      key: 'k',
      metaKey: true,
      bubbles: true
    });
    document.dispatchEvent(cmdKEvent);
    
    // Test error handling
    console.log('❌ Testing error handling...');
    document.dispatchEvent(new CustomEvent('appError', {
      detail: { message: 'Test error from modular system', type: 'warning' }
    }));
    
    console.log('✅ Modular component tests complete');
  } else {
    console.error('❌ Modular app not found');
  }
}, 2000); 