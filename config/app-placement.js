/**
 * App Placement Configuration
 * Centralized configuration for browser tab management and app placement
 */

module.exports = {
  // Single tab strategy configuration
  singleTab: {
    enabled: true,
    debugPort: 9222,
    userDataDir: '.chrome-profile', // Relative to project root
    tabTitle: 'Molecular Analysis', // why do we even want these configurations if they're already configured; don't repeat anything
    maxRetries: 3,
    
    // Chrome launch arguments
    chromeArgs: {
      // Security and permissions (updated for modern Chrome)
      security: [
        '--no-sandbox',
        '--disable-web-security',
        '--allow-running-insecure-content',
        '--use-fake-ui-for-media-stream', // Auto-grant camera permission
        '--use-fake-device-for-media-stream',
        '--disable-features=VizDisplayCompositor',
        '--ignore-certificate-errors',
        '--ignore-ssl-errors',
        '--ignore-certificate-errors-spki-list'
      ],
      
      // Performance and behavior
      performance: [
        '--disable-background-timer-throttling',
        '--disable-backgrounding-occluded-windows',
        '--disable-renderer-backgrounding',
        '--disable-background-networking',
        '--force-prefers-reduced-motion',
        '--autoplay-policy=no-user-gesture-required',
        '--disable-dev-shm-usage' // Prevents crashes on limited memory systems
      ],
      
      // UI and extensions
      ui: [
        '--no-first-run',
        '--disable-default-apps',
        '--disable-extensions'
      ]
    }
  },

  // URL strategy - HTTP first for better dev experience
  urls: {
    // Mobile devices need HTTPS for camera access only
    mobile: {
      preferHttps: true,
      requiresCamera: true,
      protocol: 'https',
      defaultPort: 3001,
      // Try HTTP first, fallback to HTTPS if camera needed
      allowHttpFallback: false
    },
    
    // Desktop always uses HTTP for better dev experience
    desktop: {
      preferHttps: false,
      requiresCamera: false,
      protocol: 'http',
      defaultPort: 8080,
      // Force HTTP even if HTTPS available
      forceHttp: true
    }
  },

  // Device detection
  deviceDetection: {
    mobileUserAgents: [
      'Mobile',
      'Android',
      'iPhone',
      'iPad',
      'iPod',
      'BlackBerry',
      'Windows Phone'
    ],
    
    // Override for testing
    forceMobile: process.env.FORCE_MOBILE === 'true',
    forceDesktop: process.env.FORCE_DESKTOP === 'true'
  },

  // Window and viewport settings
  window: {
    desktop: {
      width: 1400,
      height: 900,
      devTools: process.env.NODE_ENV === 'development'
    },
    
    mobile: {
      // Mobile uses device defaults
      devTools: false
    }
  },

  // Fallback behavior
  fallback: {
    useSystemBrowser: true,
    showManualUrl: true,
    retryAttempts: 2,
    retryDelay: 1000
  },

  // PWA and security settings
  pwa: {
    manifestPath: '/manifest.json',
    serviceWorkerPath: '/sw.js',
    requiredPermissions: ['camera'],
    
    // Camera-specific settings for mobile
    camera: {
      facingMode: 'environment',
      width: { ideal: 1920, min: 720 },
      height: { ideal: 1080, min: 480 },
      aspectRatio: { ideal: 16/9 },
      frameRate: { ideal: 30, min: 15 }
    }
  },

  // Development vs production settings
  environment: {
    development: {
      autoOpenBrowser: true,
      showDebugInfo: true,
      cleanupOnStart: true,
      reuseTabs: true
    },

    production: { // Make these not even exposed
      autoOpenBrowser: false,
      showDebugInfo: false,
      cleanupOnStart: false,
      reuseTabs: false
    }
  }
};
