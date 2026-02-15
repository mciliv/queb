// Comprehensive App Validation Rules
// These rules ensure the entire application works correctly across all test suites

const APP_VALIDATION_RULES = {
  // Core Application Requirements
  core: {
    server: {
      startup: { maxTime: 5000, healthEndpoints: ['/'] },
      apis: {
        '/analyze-text': { method: 'POST', requiredFields: ['object'], timeout: 15000 },
        '/image-molecules': { method: 'POST', requiredFields: ['imageBase64'], timeout: 20000 },
        '/generate-sdfs': { method: 'POST', requiredFields: ['smiles'], timeout: 5000 },
        '/validate-payment': { method: 'POST', requiredFields: ['device_token'], timeout: 3000 }
      },
      performance: { maxResponseTime: 2000, maxMemoryUsage: 1073741824 }
    },
    
    frontend: {
      components: {
        textInput: ['#object-input', '.validation-message'],
        camera: ['#video-feed', '#capture-btn', '.camera-status'],
        payment: ['.payment-setup', '.account-status'],
        molecular: ['.viewer-container', '.molecules-display']
      },
      performance: { maxLoadTime: 3000, maxInteractionLatency: 100 },
      accessibility: { keyboardNav: true, ariaLabels: true, focusManagement: true }
    },
    
    database: {
      optional: true,
      fallback: { fileSystem: true, gracefulDegradation: true },
      operations: { maxLatency: 1000, backupStrategy: true }
    }
  },

  // Input Validation and Data Integrity
  validation: {
    textInput: {
      minLength: 2,
      maxLength: 500,
      invalidPatterns: [
        /^(love|hate|happy|sad|angry|joy|fear|hope|dream|idea|thought|feeling|emotion)/,
        /^(running|walking|talking|thinking|sleeping|eating|drinking)$/,
        /^[a-z]{1,2}$/,
        /^[^a-z]*$/,
        /^(asdf|qwerty|test|random|nothing|something|anything|everything)$/i,
        /^(blah|meh|hmm|ugh|oof|nah|yeah|yep|nope|ok|okay)$/i,
        /^(spam|fake|scam|bot|robot|ai|gpt)$/i
      ],
      validExamples: ['water', 'ethanol', 'coffee', 'apple', 'salt', 'sugar', 'olive oil']
    },
    
    imageInput: {
      supportedFormats: ['jpeg', 'jpg', 'png', 'webp'],
      maxSize: 10485760, // 10MB
      minSize: 1024, // 1KB
      dimensions: { min: { width: 50, height: 50 }, max: { width: 4096, height: 4096 } }
    },
    
    molecularData: {
      smilesFormat: /^[A-Za-z0-9@+\-\[\]()=#$\.]+$/,
      maxSmilesLength: 1000,
      requiredFields: ['name', 'smiles'],
      optionalFields: ['formula', 'mass', 'properties']
    }
  },

  // User Experience and Workflows
  userExperience: {
    workflows: {
      textAnalysis: {
        steps: ['input-validation', 'ai-processing', 'result-display'],
        maxDuration: 20000,
        errorRecovery: true,
        progressIndicators: true
      },
      imageAnalysis: {
        steps: ['image-upload', 'crop-selection', 'ai-processing', 'result-display'],
        maxDuration: 25000,
        errorRecovery: true,
        previewSupport: true
      },
      cameraCapture: {
        steps: ['permission-request', 'camera-init', 'frame-capture', 'image-processing'],
        maxDuration: 15000,
        safariCompatibility: true,
        permissionHandling: true
      }
    },
    
    feedback: {
      loadingStates: { required: true, animation: true },
      errorMessages: { userFriendly: true, actionable: true },
      successIndicators: { clear: true, informative: true },
      progressTracking: { accurate: true, responsive: true }
    }
  },

  // Security and Privacy
  security: {
    dataHandling: {
      inputSanitization: true,
      outputEscaping: true,
      xssProtection: true,
      sqlInjectionProtection: true
    },
    
    apiSecurity: {
      rateLimit: { enabled: true, maxRequests: 100, timeWindow: 60000 },
      authentication: { tokenValidation: true, sessionManagement: true },
      encryption: { inTransit: true, sensitiveData: true }
    },
    
    privacy: {
      dataRetention: { temporary: true, userControl: true },
      analytics: { anonymized: true, optOut: true },
      cookies: { essential: true, consent: true }
    }
  },

  // Performance and Reliability
  performance: {
    benchmarks: {
      pageLoad: { maxTime: 3000, criticalResources: ['bundle.js', 'style.css'] },
      apiResponse: { maxTime: 15000, averageTime: 5000 },
      imageProcessing: { maxTime: 20000, memoryEfficient: true },
      molecularRendering: { maxTime: 5000, interactive: true }
    },
    
    reliability: {
      errorHandling: { gracefulDegradation: true, userNotification: true },
      networkResilience: { retryLogic: true, timeoutHandling: true },
      memoryManagement: { noLeaks: true, efficientCleanup: true },
      concurrency: { threadSafe: true, sessionIsolation: true }
    }
  },

  // Browser and Device Compatibility
  compatibility: {
    browsers: {
      chrome: { minVersion: 90, fullSupport: true },
      firefox: { minVersion: 88, fullSupport: true },
      safari: { minVersion: 14, cameraWorkaround: true },
      edge: { minVersion: 90, fullSupport: true }
    },
    
    devices: {
      desktop: { fullFeatures: true, keyboardShortcuts: true },
      mobile: { responsive: true, touchOptimized: true },
      tablet: { responsive: true, hybridInput: true }
    },
    
    features: {
      webgl: { required: true, fallback: false },
      camera: { required: false, gracefulFallback: true },
      fileSystem: { required: true, permissionHandling: true }
    }
  },

  // Testing and Quality Assurance
  testing: {
    coverage: {
      unit: { minimum: 80, critical: 95 },
      integration: { minimum: 70, workflows: 90 },
      e2e: { minimum: 60, criticalPaths: 100 }
    },
    
    validation: {
      inputEdgeCases: true,
      errorScenarios: true,
      performanceLimits: true,
      securityVulnerabilities: true
    },
    
    automation: {
      continuousIntegration: true,
      performanceMonitoring: true,
      regressionTesting: true,
      accessibilityChecks: true
    }
  }
};

// Validation Helper Functions
const ValidationHelpers = {
  // Check if input matches validation rules
  validateTextInput: (text) => {
    if (!text || !text.trim()) return 'Empty input not allowed';
    
    const trimmed = text.trim();
    const rules = APP_VALIDATION_RULES.validation.textInput;
    
    if (trimmed.length < rules.minLength) return `Input too short (min ${rules.minLength})`;
    if (trimmed.length > rules.maxLength) return `Input too long (max ${rules.maxLength})`;
    
    for (const pattern of rules.invalidPatterns) {
      if (pattern.test(trimmed.toLowerCase())) {
        return 'Invalid input pattern detected';
      }
    }
    
    return null; // Valid
  },
  
  // Check if API response matches expected format
  validateApiResponse: (response, endpoint) => {
    const endpointRules = APP_VALIDATION_RULES.core.server.apis[endpoint];
    if (!endpointRules) return false;
    
    // Check response structure
    if (response.error && !response.data) return true; // Error response is valid
    if (response.data && !response.error) return true; // Success response is valid
    
    return false; // Invalid response structure
  },
  
  // Check if component meets performance benchmarks
  validatePerformance: (operation, duration) => {
    const benchmarks = APP_VALIDATION_RULES.performance.benchmarks;
    const benchmark = benchmarks[operation];
    
    if (!benchmark) return true; // No benchmark defined
    
    return duration <= benchmark.maxTime;
  },
  
  // Check if security requirements are met
  validateSecurity: (context, data) => {
    const securityRules = APP_VALIDATION_RULES.security;
    
    // Basic input sanitization check
    if (context === 'input' && typeof data === 'string') {
      return !/<script|javascript:|on\w+=/i.test(data);
    }
    
    return true; // Default to secure
  },
  
  // Check if user experience requirements are met
  validateUserExperience: (workflow, elements) => {
    const workflowRules = APP_VALIDATION_RULES.userExperience.workflows[workflow];
    if (!workflowRules) return true;
    
    // Check if required elements are present
    return elements.every(element => element !== null && element !== undefined);
  }
};

module.exports = {
  APP_VALIDATION_RULES,
  ValidationHelpers
};
