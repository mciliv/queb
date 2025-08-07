#!/usr/bin/env node

const AutoTabConnector = require('./test/utils/auto-tab-connector');

async function runSeamlessTest() {
  const testCases = process.argv.slice(2);
  
  if (testCases.length === 0) {
    console.log('Usage: node seamless-test.js <test1> [test2] [test3]...');
    console.log('Example: node seamless-test.js water ethanol "black coffee"');
    return;
  }

  try {
    for (const testCase of testCases) {
      console.log(`🧪 Testing: ${testCase}`);
      
      const page = await AutoTabConnector.executeTestSilently(testCase);
      
      // Wait for analysis
      await new Promise(resolve => setTimeout(resolve, 6000));
      
      // Take screenshot
      const screenshotPath = `test/screenshots/seamless_${testCase.replace(/\s+/g, '_')}.png`;
      await page.screenshot({ path: screenshotPath, fullPage: true });
      
      console.log(`📸 Screenshot: ${screenshotPath}`);
    }
    
    console.log('✅ All tests complete - Chrome tab remains open');
    
  } catch (error) {
    console.error('❌ Error:', error.message);
  }
}

runSeamlessTest();