#!/usr/bin/env node

const puppeteer = require('puppeteer');
const path = require('path');

async function runSingleTabMolecularTests() {
  console.log('🚀 Single Tab Molecular Testing Demo');
  console.log('   This will reuse your existing Chrome tab!');
  console.log('   =====================================\n');

  let browser;
  let page;

  const testCases = [
    { name: 'Water', input: 'water' },
    { name: 'Ethanol', input: 'ethanol' },
    { name: 'Coffee', input: 'black coffee' },
    { name: 'Wine', input: 'red wine' }
  ];

  try {
    // Launch Chrome with specific settings to avoid multiple windows
    console.log('🌟 Starting Chrome with single window mode...');
    browser = await puppeteer.launch({
      headless: false,
      defaultViewport: { width: 1400, height: 900 },
      args: [
        '--no-sandbox',
        '--disable-setuid-sandbox',
        '--disable-web-security',
        '--start-maximized',
        '--no-first-run',
        '--disable-default-apps',
        '--remote-debugging-port=9223' // Different port to avoid conflicts
      ],
      slowMo: 500,
      userDataDir: './.tmp/single-tab-chrome' // Dedicated profile
    });

    // Use the first (and likely only) page
    const pages = await browser.pages();
    page = pages[0] || await browser.newPage();
    
    console.log('✅ Chrome ready! Using single tab for all tests\n');

    // Enhanced logging
    page.on('console', msg => {
      const text = msg.text();
      if (text.includes('🧪') || text.includes('API') || text.includes('molecular')) {
        console.log(`   📱 App: ${text}`);
      }
    });

    // Navigate to molecular app
    console.log('🌐 Navigating to molecular interface...');
    await page.goto('http://localhost:3001', { waitUntil: 'networkidle0' });
    await page.waitForSelector('input[type="text"]', { timeout: 10000 });
    await new Promise(resolve => setTimeout(resolve, 2000));
    
    console.log('✅ Molecular app loaded in single tab!\n');

    // Run each test in the same tab
    for (let i = 0; i < testCases.length; i++) {
      const test = testCases[i];
      console.log(`🧬 Test ${i+1}/${testCases.length}: ${test.name}`);
      console.log(`   Input: "${test.input}"`);

      // Clear and focus input
      await page.evaluate(() => {
        const input = document.querySelector('input[type="text"]');
        if (input) {
          input.value = '';
          input.focus();
        }
      });

      await new Promise(resolve => setTimeout(resolve, 500));

      // Type with visible delay
      console.log(`   ⌨️  Typing "${test.input}"...`);
      await page.type('input[type="text"]', test.input, { delay: 150 });

      console.log('   🚀 Triggering analysis...');
      await page.keyboard.press('Enter');

      console.log('   ⏳ Waiting for molecular visualization...');
      await new Promise(resolve => setTimeout(resolve, 8000));

      // Check page content
      const hasContent = await page.evaluate(() => {
        const content = document.body.innerText;
        return {
          analyzing: content.includes('Analyzing'),
          molecular: content.includes('molecular') || content.includes('Molecular'),
          length: content.length
        };
      });

      console.log(`   📊 Content check: ${JSON.stringify(hasContent)}`);
      console.log(`   ✅ ${test.name} test complete\n`);

      // Pause for observation
      if (i < testCases.length - 1) {
        console.log('   ⏸️  Pausing 3 seconds for visual inspection...');
        await new Promise(resolve => setTimeout(resolve, 3000));
      }
    }

    console.log('🎯 All tests complete! Chrome tab will remain open for inspection.');
    console.log('   Close the browser manually when done.\n');

    // Don't close browser automatically
    // await browser.close();

  } catch (error) {
    console.error('❌ Error during testing:', error.message);
    if (browser) {
      await browser.close();
    }
  }
}

// Run the demo
runSingleTabMolecularTests().catch(console.error);