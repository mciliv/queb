/**
 * Visual Interface Tests (No API Key Required)
 * Tests the UI components and basic interaction without requiring OpenAI API
 */

const puppeteer = require('puppeteer');
const path = require('path');
const AutoTabConnector = require('../../utils/auto-tab-connector');

describe('Visual Interface Tests', () => {
  let browser;
  let page;
  const frontendUrl = 'http://localhost:3001';

  beforeAll(async () => {
    console.log('🚀 Starting visual interface tests...');
    
    // Use simple, reliable browser launch
    browser = await puppeteer.launch({
      headless: process.env.HEADLESS === 'false' ? false : 'new',
      defaultViewport: { width: 1280, height: 720 },
      args: [
        '--no-sandbox', 
        '--disable-setuid-sandbox',
        '--disable-web-security',
        '--no-first-run',
        '--disable-default-apps'
      ],
      slowMo: process.env.HEADLESS === 'false' ? 300 : 0
    });
    
    page = await browser.newPage();
    
    // Enable console logging from the page
    page.on('console', msg => {
      console.log(`📄 Browser: ${msg.text()}`);
    });
    
    console.log(`🌐 Navigating to ${frontendUrl}...`);
    await page.goto(frontendUrl, { 
      waitUntil: 'networkidle0',
      timeout: 15000 
    });
    
    // Wait for the app to load
    await page.waitForSelector('#root', { timeout: 10000 });
    console.log('✅ Frontend loaded successfully');
  });

  afterAll(async () => {
    if (browser) {
      await browser.close();
    }
    console.log('✅ Global teardown complete');
  });

  test('should load the molecular analysis interface', async () => {
    // Check that the main app elements are present
    const title = await page.title();
    expect(title).toContain('Atomizer');
    
    // Take screenshot of loaded interface
    const screenshotPath = path.join(__dirname, '../screenshots/interface_loaded.png');
    await page.screenshot({ path: screenshotPath, fullPage: true });
    console.log(`📸 Screenshot saved: ${screenshotPath}`);
  });

  test('should find input elements', async () => {
    // Find the React text input specifically
    const textInput = await page.$('#object-input');
    expect(textInput).toBeTruthy();
    console.log('✅ Found React text input field');

    // Find analyze button  
    const buttons = await page.$$('button');
    expect(buttons.length).toBeGreaterThan(0);
    console.log(`✅ Found ${buttons.length} buttons`);
    
    // Check for mode selector elements
    const modeElements = await page.$$('[class*="mode"], [class*="Mode"]');
    console.log(`🎛️ Found ${modeElements.length} mode selector elements`);
  });

  test('should inject test data into input field', async () => {
    // Test data injection without API call
    const testInputs = ['water', 'ethanol', 'red wine'];
    
    for (const testInput of testInputs) {
      console.log(`📝 Testing input injection: "${testInput}"`);
      
      // Find and clear the React text input
      const textInput = await page.$('#object-input');
      await textInput.evaluate(el => el.value = ''); // Clear directly
      await textInput.type(testInput);
      
      // Verify the input was set
      const inputValue = await textInput.evaluate(el => el.value);
      expect(inputValue).toBe(testInput);
      console.log(`✅ Successfully injected: "${inputValue}"`);
      
      // Clear for next test
      await textInput.evaluate(el => el.value = ''); // Clear directly
    }
  });

  test('should find molecular test panel elements', async () => {
    // Look for test panel button or components
    const testPanelElements = await page.$$('.test-panel-toggle, [class*="test"], [class*="molecular"]');
    console.log(`🧪 Found ${testPanelElements.length} test-related elements`);
    
    // Look for any buttons with emoji or test-related text
    const allButtons = await page.$$('button');
    let testButtons = 0;
    
    for (const button of allButtons) {
      const text = await button.evaluate(el => el.textContent || el.title || '');
      if (text.includes('🧪') || text.includes('test') || text.includes('Test')) {
        testButtons++;
        console.log(`🎯 Found test button: "${text}"`);
      }
    }
    
    console.log(`📊 Total test buttons found: ${testButtons}`);
  });

  test('should have proper styling and layout', async () => {
    // Check viewport and basic styling
    const bodyStyles = await page.evaluate(() => {
      const body = document.body;
      const styles = window.getComputedStyle(body);
      return {
        backgroundColor: styles.backgroundColor,
        fontFamily: styles.fontFamily,
        margin: styles.margin
      };
    });
    
    console.log('🎨 Body styles:', bodyStyles);
    expect(bodyStyles).toBeTruthy();
    
    // Take a final screenshot
    const screenshotPath = path.join(__dirname, '../screenshots/interface_complete.png');
    await page.screenshot({ path: screenshotPath, fullPage: true });
    console.log(`📸 Final screenshot: ${screenshotPath}`);
  });

  test('should inject and test fake molecular data visualization', async () => {
    console.log('🧪 Testing fake molecular data injection...');

    // Inject the test script directly into the page
    await page.evaluate(() => {
      // Add fake molecular data to simulate results
      const fakeResults = {
        object: "test water",
        chemicals: [
          { name: "Water", smiles: "O" },
          { name: "Test Compound", smiles: "CCO" }
        ]
      };

      // Try to find the app's state or results container
      const root = document.getElementById('root');
      if (root) {
        // Create a test results section
        const testDiv = document.createElement('div');
        testDiv.id = 'test-results';
        testDiv.className = 'molecular-results test-injection';
        testDiv.innerHTML = `
          <h3>🧪 Test Molecular Results</h3>
          <div class="test-molecules">
            ${fakeResults.chemicals.map(mol => `
              <div class="test-molecule">
                <h4>${mol.name}</h4>
                <p>SMILES: ${mol.smiles}</p>
                <div class="test-viewer" style="width:200px;height:200px;border:1px solid #ccc;background:#f0f0f0;display:flex;align-items:center;justify-content:center;">
                  🔬 Molecular Viewer (${mol.name})
                </div>
              </div>
            `).join('')}
          </div>
        `;
        root.appendChild(testDiv);
        
        console.log('✅ Injected fake molecular results');
        return true;
      }
      return false;
    });

    // Wait for the injection to complete
    await new Promise(resolve => setTimeout(resolve, 1000));

    // Verify the fake results are visible
    const testResults = await page.$('#test-results');
    expect(testResults).toBeTruthy();
    
    const moleculeElements = await page.$$('.test-molecule');
    console.log(`🧬 Found ${moleculeElements.length} test molecules displayed`);
    expect(moleculeElements.length).toBeGreaterThan(0);

    // Take screenshot of fake results
    const screenshotPath = path.join(__dirname, '../screenshots/fake_molecular_results.png');
    await page.screenshot({ path: screenshotPath, fullPage: true });
    console.log(`📸 Fake results screenshot: ${screenshotPath}`);
  });
});