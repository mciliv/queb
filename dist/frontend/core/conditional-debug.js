// Conditional debug loader - only loads debug utilities in development
(function() {
  const IS_PRODUCTION = window.location.hostname !== 'localhost' && 
                        window.location.hostname !== '127.0.0.1' &&
                        !window.location.hostname.includes('dev') &&
                        !window.location.search.includes('debug=true');

  // Set global production flag
  window.IS_PRODUCTION = IS_PRODUCTION;

  
    };
  } else {
    // Production mode - create stub functions
    window.debug = {
      log: function() {},
      warn: function() {},
      error: function() {},
      inspect: function() {},
      highlightElement: function() {},
      showBoundingBoxes: function() {},
      hideBoundingBoxes: function() {}
    };
    
    // Override console methods in production
    if (window.location.hostname !== 'localhost') {
      const noop = function() {};
      console.log = noop;
      console.debug = noop;
      console.info = noop;
      // Keep warn and error for critical issues
    }
  }
})(); 