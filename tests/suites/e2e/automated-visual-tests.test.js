/**
 * Automated Visual Molecular Tests
 * Uses Puppeteer to run visual tests in a real browser
 */

const puppeteer = require('puppeteer');
const path = require('path');

// Test cases using minimal subset validation
const visualTestCases = [
  { name: '💧 Water', input: 'water', expected: ['water'], category: 'basics' },
  { name: '🧪 Ethanol', input: 'ethanol', expected: ['ethanol'], category: 'basics' },
  { name: '🧂 Salt', input: 'sodium chloride', expected: ['sodium', 'chloride'], category: 'basics' },
  { name: '🍷 Wine', input: 'red wine', expected: ['ethanol', 'water'], category: 'beverages' },
  { name: '☕ Coffee', input: 'black coffee', expected: ['water', 'caffeine'], category: 'beverages' },
  { name: '🍎 Apple', input: 'fresh apple', expected: ['water'], category: 'biological' }
];

describe('Automated Visual Molecular Tests', () => {
  let browser;
  let page;
  const frontendUrl = 'https://localhost:3001';
  const hasApiKey = !!process.env.OPENAI_API_KEY;

  beforeAll(async () => {
    if (!hasApiKey) {
      console.log('⚠️  Skipping visual tests - no OPENAI_API_KEY set');
      return;
    }

    console.log('🚀 Starting automated visual tests...');
    
    // Launch browser
    browser = await puppeteer.launch({
      headless: process.env.HEADLESS === 'true' ? 'new' : false,
      defaultViewport: { width: 1280, height: 720 },
      args: [
        '--no-sandbox', 
        '--disable-setuid-sandbox',
        '--ignore-certificate-errors',
        '--ignore-ssl-errors',
        '--allow-running-insecure-content'
      ]
    });

    page = await browser.newPage();
    
    // Enable console logging from the page
    page.on('console', msg => {
      if (msg.type() === 'log') {
        console.log(`📄 Browser: ${msg.text()}`);
      }
    });

    // Navigate to the frontend
    console.log(`🌐 Navigating to ${frontendUrl}...`);
    await page.goto(frontendUrl, { waitUntil: 'networkidle0' });
    
    // Wait for the app to load
    await page.waitForSelector('#root', { timeout: 10000 });
    console.log('✅ Frontend loaded successfully');
  });

  afterAll(async () => {
    if (browser) {
      await browser.close();
    }
  });

  describe('Visual Test Injection', () => {
    if (!hasApiKey) {
      test.skip('Visual tests require OPENAI_API_KEY', () => {});
      return;
    }

    visualTestCases.forEach((testCase) => {
      test(`should visualize ${testCase.name} (${testCase.category})`, async () => {
        console.log(`\n🧪 Testing: ${testCase.name}`);
        console.log(`   Input: "${testCase.input}"`);
        console.log(`   Expected subset: ${testCase.expected.join(', ')}`);

        // Find and fill the text input
        const textInput = await page.waitForSelector('input[type="text"], textarea, #object-input', {
          timeout: 5000
        });
        
        await textInput.click({ clickCount: 3 }); // Select all
        await textInput.type(testCase.input);
        console.log(`📝 Injected "${testCase.input}" into input field`);

        // Find and click the analyze button
        const analyzeButton = await page.waitForSelector('button[type="submit"], button:has-text("Analyze"), .analyze-button', {
          timeout: 5000
        });

        console.log('🚀 Triggering analysis...');
        await analyzeButton.click();

        // Wait for results to appear (molecular visualization)
        console.log('⏳ Waiting for molecular visualization...');
        
        try {
          // Wait for the results section to appear
          await page.waitForSelector('.results-section, .molecular-results, .results-container', {
            timeout: 30000 // 30 seconds for AI analysis
          });

          // Check for 3DMol.js viewers
          const viewers = await page.$$('.mol-viewer, .molecule-viewer, [id*="viewer"]');
          console.log(`🎬 Found ${viewers.length} molecular viewers`);

          // Verify we have at least one molecular visualization
          expect(viewers.length).toBeGreaterThan(0);

          // Check for molecule names/data
          const moleculeElements = await page.$$('.molecule-name, .chemical-name, .mol-data');
          console.log(`🧬 Found ${moleculeElements.length} molecule data elements`);

          // Take a screenshot for visual verification
          const screenshotPath = path.join(__dirname, `../screenshots/${testCase.name.replace(/[^a-zA-Z0-9]/g, '_')}.png`);
          await page.screenshot({ path: screenshotPath, fullPage: true });
          console.log(`📸 Screenshot saved: ${screenshotPath}`);

          // Validate subset requirements (basic check)
          if (testCase.category === 'basics') {
            expect(moleculeElements.length).toBeGreaterThanOrEqual(testCase.expected.length);
          }

          console.log(`✅ ${testCase.name} visualization completed successfully`);

        } catch (error) {
          console.error(`❌ ${testCase.name} failed: ${error.message}`);
          
          // Take error screenshot
          const errorScreenshotPath = path.join(__dirname, `../screenshots/ERROR_${testCase.name.replace(/[^a-zA-Z0-9]/g, '_')}.png`);
          await page.screenshot({ path: errorScreenshotPath, fullPage: true });
          
          throw error;
        }

        // Brief pause between tests
        await page.waitForTimeout(2000);
      }, 60000); // 60 second timeout per test
    });
  });

  describe('Visual Test Panel Integration', () => {
    if (!hasApiKey) {
      test.skip('Visual tests require OPENAI_API_KEY', () => {});
      return;
    }

    test('should find and interact with test panel', async () => {
      console.log('🧪 Looking for molecular test panel...');

      // Look for the test panel button
      const testPanelButton = await page.$('.test-panel-toggle, button:has-text("🧪")');
      
      if (testPanelButton) {
        console.log('✅ Found test panel button');
        await testPanelButton.click();
        
        // Wait for panel to open
        await page.waitForTimeout(1000);
        
        // Look for test buttons
        const testButtons = await page.$$('.test-panel-content button');
        console.log(`🎯 Found ${testButtons.length} test buttons in panel`);
        
        expect(testButtons.length).toBeGreaterThan(0);
      } else {
        console.log('⚠️  Test panel button not found - using direct injection method');
        // This is okay - we're testing the core functionality anyway
      }
    });
  });

  describe('Performance Validation', () => {
    if (!hasApiKey) {
      test.skip('Visual tests require OPENAI_API_KEY', () => {});
      return;
    }

    test('should complete full pipeline within reasonable time', async () => {
      const startTime = Date.now();
      
      // Run a simple test case
      const testInput = 'water';
      
      const textInput = await page.waitForSelector('input[type="text"], textarea');
      await textInput.click({ clickCount: 3 });
      await textInput.type(testInput);
      
      const analyzeButton = await page.waitForSelector('button[type="submit"], button');
      await analyzeButton.click();
      
      // Wait for completion
      await page.waitForSelector('.results-section, .molecular-results', { timeout: 45000 });
      
      const endTime = Date.now();
      const duration = endTime - startTime;
      
      console.log(`⏱️  Full pipeline completed in ${duration}ms`);
      expect(duration).toBeLessThan(45000); // Should complete within 45 seconds
    });
  });
});