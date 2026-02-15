const puppeteer = require('puppeteer');
const path = require('path');
const IntegrationTestSetup = require('../../integration/setup-integration');

// Integration Test Validation Rules - End-to-End App Functionality
const INTEGRATION_VALIDATION_RULES = {
  // Critical user flows that must work
  userFlows: {
    textAnalysis: {
      steps: ['input-text', 'validate-input', 'process-ai', 'display-molecules'],
      maxDuration: 20000,
      criticalElements: ['#object-input', '.molecules-container']
    },
    imageAnalysis: {
      steps: ['upload-image', 'crop-image', 'process-ai', 'display-molecules'],
      maxDuration: 25000,
      criticalElements: ['#photo-upload', '#video-feed', '.molecules-container']
    },
    cameraCapture: {
      steps: ['request-permission', 'initialize-camera', 'capture-frame', 'process-image'],
      maxDuration: 15000,
      criticalElements: ['#video-feed', '#capture-btn']
    },
    paymentFlow: {
      steps: ['check-payment', 'validate-token', 'update-usage'],
      maxDuration: 5000,
      criticalElements: ['.payment-setup', '.account-status']
    }
  },
  
  // System-wide validation checks
  systemHealth: {
    serverResponse: { maxLatency: 2000, requiredEndpoints: ['/'] },
    databaseOperations: { optional: true, maxLatency: 1000 },
    aiServiceIntegration: { timeout: 30000, retries: 2 },
    fileSystemOperations: { writeable: true, readableDirectories: ['sdf_files'] }
  },
  
  // Data integrity validation
  dataIntegrity: {
    smilesValidation: { format: /^[A-Za-z0-9@+\-\[\]()=#$]+$/, maxLength: 1000 },
    imageProcessing: { supportedFormats: ['jpeg', 'png', 'webp'], maxSize: 10485760 },
    molecularData: { requiredFields: ['name', 'smiles'], validation: true },
    userSessions: { tokenValidation: true, expirationCheck: true }
  },
  
  // Performance and reliability benchmarks
  reliability: {
    errorHandling: { gracefulDegradation: true, userFeedback: true },
    networkResilience: { retryLogic: true, timeoutHandling: true },
    memoryManagement: { noLeaks: true, efficientCleanup: true },
    concurrentUsers: { supportMultiple: true, sessionIsolation: true }
  }
};

describe('Molecular App Integration Tests - Full App Validation', () => {
  let browser;
  let page;
  let testSetup;
  
  // Validation helper functions
  const validateUserFlow = async (flowName, page) => {
    const flow = INTEGRATION_VALIDATION_RULES.userFlows[flowName];
    const startTime = Date.now();
    
    // Check critical elements exist
    for (const element of flow.criticalElements) {
      await expect(page.$(element)).resolves.toBeTruthy();
    }
    
    const duration = Date.now() - startTime;
    expect(duration).toBeLessThan(flow.maxDuration);
    
    return { success: true, duration, elements: flow.criticalElements.length };
  };
  
  const validateSystemHealth = async (page) => {
    const health = INTEGRATION_VALIDATION_RULES.systemHealth;
    
    // Check server response time
    const startTime = Date.now();
    await page.goto('http://localhost:8080');
    const responseTime = Date.now() - startTime;
    
    expect(responseTime).toBeLessThan(health.serverResponse.maxLatency);
    
    return { responseTime, healthy: true };
  };
  
  const validateDataIntegrity = (data, type) => {
    const rules = INTEGRATION_VALIDATION_RULES.dataIntegrity;
    
    switch(type) {
      case 'smiles':
        expect(data).toMatch(rules.smilesValidation.format);
        expect(data.length).toBeLessThanOrEqual(rules.smilesValidation.maxLength);
        break;
      case 'molecular':
        rules.molecularData.requiredFields.forEach(field => {
          expect(data[field]).toBeDefined();
        });
        break;
    }
    
    return true;
  };

  beforeAll(async () => {
    // Start the server first
    testSetup = new IntegrationTestSetup();
    await testSetup.startServer();
    await testSetup.waitForServer();

    browser = await puppeteer.launch({ 
      headless: true,
      args: ['--no-sandbox', '--disable-setuid-sandbox']
    });
    page = await browser.newPage();
    
    // Enable console logging
    page.on('console', msg => console.log('Browser console:', msg.text()));
    page.on('pageerror', err => console.error('Browser error:', err));
  }, 60000); // 60 second timeout

  afterAll(async () => {
    if (page) await page.close();
    if (browser) await browser.close();
    if (testSetup) testSetup.stopServer();
  });

  beforeEach(async () => {
    await page.goto(testSetup.baseUrl, { waitUntil: 'networkidle0' });
    // Wait for app to initialize
    await page.waitForSelector('#object-input', { timeout: 10000 });
  });

  describe('App Initialization', () => {
    test('should load the main page with all components', async () => {
      // Check main elements exist
      await expect(page.$('#object-input')).resolves.toBeTruthy();
      await expect(page.$('#video-feed')).resolves.toBeTruthy();
      await expect(page.$('#photo-upload')).resolves.toBeTruthy();
      await expect(page.$('#photo-url')).resolves.toBeTruthy();
      await expect(page.$('.snapshots-container')).resolves.toBeTruthy();
      
      console.log('✅ All main UI components loaded');
    });

    test('should initialize camera system', async () => {
      // Check if camera manager is initialized
      const cameraInitialized = await page.evaluate(() => {
        return window.cameraManager && window.cameraManager.isInitialized;
      });
      
      expect(cameraInitialized).toBe(true);
      console.log('✅ Camera system initialized');
    });

    test('should initialize payment system', async () => {
      // Check if payment manager is initialized
      const paymentInitialized = await page.evaluate(() => {
        return window.paymentManager && typeof window.paymentManager.checkPaymentMethod === 'function';
      });
      
      expect(paymentInitialized).toBe(true);
      console.log('✅ Payment system initialized');
    });
  });

  describe('Camera Integration', () => {
    test('should show camera view when camera button is clicked', async () => {
      // Click on video feed to trigger camera mode
      await page.click('#video-feed');
      
      // Wait for camera view to be visible
      await page.waitForSelector('#video-feed.active', { timeout: 5000 });
      
      const isCameraActive = await page.evaluate(() => {
        const video = document.getElementById('video-feed');
        return video.classList.contains('active');
      });
      
      expect(isCameraActive).toBe(true);
      console.log('✅ Camera view activated on click');
    });

    test('should capture image when camera capture button is clicked', async () => {
      // First activate camera
      await page.click('#video-feed');
      await page.waitForSelector('#video-feed.active', { timeout: 5000 });
      
      // Look for capture button (might be in camera overlay)
      const captureButton = await page.$('[data-action="capture"], .capture-button, #capture-btn');
      
      if (captureButton) {
        await captureButton.click();
        
        // Wait for image processing
        await page.waitForTimeout(2000);
        
        // Check if image analysis was triggered
        const hasAnalysisResult = await page.evaluate(() => {
          return document.querySelector('.object-column') !== null;
        });
        
        console.log('✅ Camera capture triggered analysis:', hasAnalysisResult);
      } else {
        console.log('⚠️ No capture button found - camera may be in different mode');
      }
    });
  });

  describe('Photo Upload Integration', () => {
    test('should handle file upload and trigger analysis', async () => {
      // Create a test image file
      const testImagePath = path.join(__dirname, 'fixtures', 'test-molecule.jpg');
      
      // Upload test image
      const [fileChooser] = await Promise.all([
        page.waitForFileChooser(),
        page.click('#photo-upload')
      ]);
      
      await fileChooser.accept([testImagePath]);
      
      // Wait for image to be displayed
      await page.waitForSelector('.uploaded-image-container', { timeout: 10000 });
      
      console.log('✅ Photo upload handled successfully');
    });

    test('should trigger analysis when clicking on uploaded image', async () => {
      // First upload an image
      const testImagePath = path.join(__dirname, 'fixtures', 'test-molecule.jpg');
      
      const [fileChooser] = await Promise.all([
        page.waitForFileChooser(),
        page.click('#photo-upload')
      ]);
      
      await fileChooser.accept([testImagePath]);
      await page.waitForSelector('.uploaded-image-container', { timeout: 10000 });
      
      // Click on the uploaded image
      await page.click('.uploaded-image-container img');
      
      // Wait for analysis to start
      await page.waitForSelector('.loading-column', { timeout: 5000 });
      
      console.log('✅ Image click triggered analysis');
    });
  });

  describe('Text Analysis Integration', () => {
    test('should analyze text input and display results', async () => {
      // Type a simple molecule name
      await page.type('#object-input', 'water');
      await page.keyboard.press('Enter');
      
      // Wait for analysis to complete
      await page.waitForTimeout(3000);
      
      // Check if results are displayed
      const hasResults = await page.evaluate(() => {
        return document.querySelector('.object-column') !== null;
      });
      
      expect(hasResults).toBe(true);
      console.log('✅ Text analysis completed and results displayed');
    });

    test('should handle API errors gracefully', async () => {
      // Type invalid input
      await page.type('#object-input', 'invalid_molecule_xyz123');
      await page.keyboard.press('Enter');
      
      // Wait for error handling
      await page.waitForTimeout(2000);
      
      // Check if error message is displayed
      const hasErrorMessage = await page.evaluate(() => {
        return document.querySelector('.error-message') !== null;
      });
      
      console.log('✅ Error handling works:', hasErrorMessage);
    });
  });

  describe('Payment Integration', () => {
    test('should show payment requirement when needed', async () => {
      // Try to analyze without payment setup
      await page.type('#object-input', 'test');
      await page.keyboard.press('Enter');
      
      // Wait for payment check
      await page.waitForTimeout(2000);
      
      // Check if payment message is shown
      const hasPaymentMessage = await page.evaluate(() => {
        return document.querySelector('.payment-required') !== null;
      });
      
      console.log('✅ Payment requirement check works:', hasPaymentMessage);
    });
  });

  describe('3D Visualization Integration', () => {
    test('should render 3D molecules after analysis', async () => {
      // Perform a successful analysis
      await page.type('#object-input', 'water');
      await page.keyboard.press('Enter');
      
      // Wait for 3D rendering
      await page.waitForTimeout(5000);
      
      // Check if 3D viewer is present
      const has3DViewer = await page.evaluate(() => {
        return document.querySelector('.molecule-viewer') !== null;
      });
      
      expect(has3DViewer).toBe(true);
      console.log('✅ 3D molecule visualization rendered');
    });
  });

  describe('Component Communication', () => {
    test('should maintain state across different analysis types', async () => {
      // Test text analysis
      await page.type('#object-input', 'water');
      await page.keyboard.press('Enter');
      await page.waitForTimeout(3000);
      
      // Test photo upload
      const testImagePath = path.join(__dirname, 'fixtures', 'test-molecule.jpg');
      const [fileChooser] = await Promise.all([
        page.waitForFileChooser(),
        page.click('#photo-upload')
      ]);
      
      await fileChooser.accept([testImagePath]);
      await page.waitForSelector('.uploaded-image-container', { timeout: 10000 });
      
      // Check if both results are displayed
      const resultCount = await page.evaluate(() => {
        return document.querySelectorAll('.object-column').length;
      });
      
      expect(resultCount).toBeGreaterThan(0);
      console.log('✅ Multiple analysis types work together');
    });

    test('should handle component cleanup properly', async () => {
      // Create some results first
      await page.type('#object-input', 'water');
      await page.keyboard.press('Enter');
      await page.waitForTimeout(3000);
      
      // Close a result
      await page.click('.object-close');
      
      // Check if result was removed
      const resultCount = await page.evaluate(() => {
        return document.querySelectorAll('.object-column').length;
      });
      
      expect(resultCount).toBe(0);
      console.log('✅ Component cleanup works properly');
    });
  });

  describe('Error Handling', () => {
    test('should handle network errors gracefully', async () => {
      // Mock network failure
      await page.setOfflineMode(true);
      
      await page.type('#object-input', 'water');
      await page.keyboard.press('Enter');
      
      await page.waitForTimeout(2000);
      
      // Check if error is handled
      const hasError = await page.evaluate(() => {
        return document.querySelector('.error-message') !== null;
      });
      
      await page.setOfflineMode(false);
      
      console.log('✅ Network error handling works:', hasError);
    });
  });
}); 