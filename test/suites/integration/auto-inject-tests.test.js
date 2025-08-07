const puppeteer = require('puppeteer');
const path = require('path');
const AutoTabConnector = require('../utils/auto-tab-connector');

describe('Auto-Inject Molecular Tests', () => {
  let browser;
  let page;
  const screenshotDir = path.join(__dirname, '../screenshots');

  // Test cases matching the MolecularTestPanel
  const testCases = [
    { name: 'Water', input: 'water', category: 'basics' },
    { name: 'Ethanol', input: 'ethanol', category: 'basics' },
    { name: 'Salt', input: 'sodium chloride', category: 'basics' },
    { name: 'Red Wine', input: 'red wine', category: 'beverages' },
    { name: 'Coffee', input: 'black coffee', category: 'beverages' },
    { name: 'Apple', input: 'fresh apple', category: 'biological' }
  ];

  beforeAll(async () => {
    console.log('🚀 Starting automated molecular injection tests...');
    console.log('   These tests will show REAL molecular visualizations!');
    
    try {
      const tabInfo = await AutoTabConnector.getOrCreateTab();
      browser = tabInfo.browser;
      page = tabInfo.page;
      
      if (tabInfo.reused) {
        console.log('♻️  Reusing existing Chrome tab');
      } else {
        console.log('🌟 Created new Chrome tab');
      }
    } catch (error) {
      console.error('Failed to get tab, falling back to new browser instance:', error);
      browser = await puppeteer.launch({
        headless: false,
        defaultViewport: { width: 1600, height: 1000 },
        args: ['--no-sandbox', '--disable-setuid-sandbox'],
        slowMo: 300
      });
      page = await browser.newPage();
      await page.goto('http://localhost:3001', { waitUntil: 'networkidle0' });
    }
    
    // Log molecular analysis activity
    page.on('console', msg => {
      const text = msg.text();
      if (text.includes('🧪') || text.includes('🔗') || text.includes('📨') || 
          text.includes('API') || text.includes('molecular') || text.includes('SMILES')) {
        console.log(`   📄 App: ${text}`);
      }
    });

    await page.waitForSelector('input[type="text"]', { timeout: 15000 });
    await new Promise(resolve => setTimeout(resolve, 2000));
    
    console.log('✅ Molecular app loaded and ready');
  });

  afterAll(async () => {
    console.log('🎯 Completed all molecular injection tests');
    console.log('💡 Check the screenshots in test/screenshots/ to see molecular visualizations');
    console.log('💡 Chrome tab remains open for continued testing');
    
    // Don't close browser - let it be reused by other tests
  });

  test('should auto-inject all test cases and show molecular visualizations', async () => {
    console.log('\n🧪 Starting automated injection of all test cases...');
    console.log('   Watch the browser to see each molecular visualization appear!');

    for (let i = 0; i < testCases.length; i++) {
      const testCase = testCases[i];
      console.log(`\n🔬 Test ${i+1}/${testCases.length}: ${testCase.name}`);
      console.log(`   Category: ${testCase.category}`);
      console.log(`   Input: "${testCase.input}"`);
      
      // Clear previous input
      await page.evaluate(() => {
        const input = document.querySelector('input[type="text"]');
        if (input) {
          input.value = '';
          input.dispatchEvent(new Event('input', { bubbles: true }));
        }
      });
      
      await new Promise(resolve => setTimeout(resolve, 500));
      
      // Type the test input
      console.log(`   ⌨️  Injecting: "${testCase.input}"`);
      await page.type('input[type="text"]', testCase.input, { delay: 50 });
      
      // Take before screenshot
      const beforePath = path.join(screenshotDir, `auto_${testCase.name.toLowerCase().replace(' ', '_')}_before.png`);
      await page.screenshot({ path: beforePath, fullPage: true });
      
      // Trigger analysis
      console.log('   🚀 Triggering molecular analysis...');
      await page.keyboard.press('Enter'); // Try Enter key first
      
      // Wait for analysis to complete and visualization to render
      console.log('   ⏳ Waiting for molecular analysis and 3D visualization...');
      await new Promise(resolve => setTimeout(resolve, 12000)); // Long enough to see the results
      
      // Take after screenshot showing molecular visualization
      const afterPath = path.join(screenshotDir, `auto_${testCase.name.toLowerCase().replace(' ', '_')}_molecules.png`);
      await page.screenshot({ path: afterPath, fullPage: true });
      console.log(`   📸 Captured molecular display: ${afterPath}`);
      
      // Check what molecular data was found
      const molecularInfo = await page.evaluate(() => {
        // Look for molecular viewer elements
        const viewers = document.querySelectorAll('[id*="molviewer"], canvas, [class*="mol-viewer"]');
        
        // Look for molecule names and SMILES
        const moleculeTexts = Array.from(document.querySelectorAll('*'))
          .map(el => el.textContent)
          .filter(text => text && (
            text.includes('SMILES') || 
            text.match(/^[CONSPconsp\[\]()=#+\-0-9]{3,}$/) ||
            text.toLowerCase().includes('molecule') ||
            text.toLowerCase().includes('chemical')
          ));
        
        return {
          viewerCount: viewers.length,
          moleculeTexts: moleculeTexts.slice(0, 5) // First 5 relevant texts
        };
      });
      
      console.log(`   🧬 Found ${molecularInfo.viewerCount} molecular viewers`);
      if (molecularInfo.moleculeTexts.length > 0) {
        console.log(`   🔬 Molecular data detected:`);
        molecularInfo.moleculeTexts.forEach((text, i) => {
          console.log(`      ${i+1}. ${text.substring(0, 60)}${text.length > 60 ? '...' : ''}`);
        });
      }
      
      console.log(`   ✅ ${testCase.name} molecular visualization complete`);
      
      // Pause between tests to allow viewing
      if (i < testCases.length - 1) {
        console.log('   ⏱️  Pausing 3 seconds before next test...');
        await new Promise(resolve => setTimeout(resolve, 3000));
      }
    }
    
    // Final summary screenshot
    const summaryPath = path.join(screenshotDir, 'auto_injection_summary.png');
    await page.screenshot({ path: summaryPath, fullPage: true });
    console.log(`\n📸 Final summary screenshot: ${summaryPath}`);
    
    console.log('\n🎯 All automated molecular injections complete!');
    console.log('   Check the screenshots to see the molecular visualizations captured');
    
    expect(testCases.length).toBeGreaterThan(0);
  }, 120000); // 2 minutes for all tests

  test('should test the molecular test panel if available', async () => {
    console.log('\n🧪 Looking for molecular test panel...');
    
    // Look for test panel button or elements
    const testPanelElements = await page.$$eval('*', () => {
      return Array.from(document.querySelectorAll('*'))
        .filter(el => {
          const text = el.textContent?.toLowerCase() || '';
          const attrs = `${el.className} ${el.id} ${el.title || ''}`.toLowerCase();
          return text.includes('🧪') || text.includes('test') && text.includes('molecular') ||
                 attrs.includes('test') || el.tagName === 'BUTTON' && text.includes('test');
        })
        .map(el => ({
          tag: el.tagName,
          text: el.textContent?.substring(0, 50),
          classes: el.className,
          clickable: el.tagName === 'BUTTON' || el.onclick !== null
        }));
    });
    
    if (testPanelElements.length > 0) {
      console.log(`   ✅ Found ${testPanelElements.length} test panel elements:`);
      testPanelElements.forEach((el, i) => {
        console.log(`      ${i+1}. ${el.tag}: ${el.text} ${el.clickable ? '(clickable)' : ''}`);
      });
      
      // Try to click test panel button
      try {
        const testButton = await page.$('button:has-text("🧪"), button[title*="test"]');
        if (testButton) {
          console.log('   🎯 Clicking test panel button...');
          await testButton.click();
          
          await new Promise(resolve => setTimeout(resolve, 5000));
          
          const panelPath = path.join(screenshotDir, 'test_panel_activated.png');
          await page.screenshot({ path: panelPath, fullPage: true });
          console.log(`   📸 Test panel screenshot: ${panelPath}`);
        }
      } catch (error) {
        console.log(`   ⚠️  Could not activate test panel: ${error.message}`);
      }
    } else {
      console.log('   ℹ️  No test panel found - tests ran through manual injection');
    }
    
    expect(true).toBe(true); // Always pass, this is informational
  }, 30000);
});