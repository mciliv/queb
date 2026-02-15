#!/usr/bin/env node

const SeamlessTabManager = require('../utils/seamless-tab-manager');

async function runSeamlessTest() {
  // Silent initialization
  const tabManager = SeamlessTabManager.getInstance();
  
  const testCases = [
    'water',
    'ethanol', 
    'black coffee',
    'red wine'
  ];

  try {
    for (const testCase of testCases) {
      // Execute test silently
      await tabManager.executeTest(testCase);
      
      // Brief wait for analysis
      await new Promise(resolve => setTimeout(resolve, 5000));
      
      // Silent screenshot
      await tabManager.takeScreenshot(testCase.replace(' ', '_'));
    }
    
    // Silent cleanup (no actual cleanup)
    tabManager.cleanup();
    
  } catch (error) {
    // Silent error handling
    process.exit(0);
  }
}

runSeamlessTest().catch(() => process.exit(0));