#!/usr/bin/env node
/**
 * Test Injection Script - Triggers analysis in live running app
 * 
 * This script makes HTTP requests to a running Queb server to trigger
 * molecular analysis, simulating user input but via API calls.
 * Perfect for development workflows where you want to see results
 * immediately in the live app.
 */

const fetch = require('node-fetch').default || require('node-fetch');

const SERVER_URL = 'http://localhost:8080';

// Test case definitions - what to analyze
const TEST_CASES = {
  coffee: {
    input: "coffee",
    description: "☕ Coffee with caffeine, water, chlorogenic acid",
    lookupMode: "GPT-5"
  },

  saltwater: {
    input: "salt water",
    description: "🧂 Simple ionic solution - sodium chloride + water",
    lookupMode: "GPT-5"
  },

  aspirin: {
    input: "aspirin tablet",
    description: "💊 Common medication - acetylsalicylic acid",
    lookupMode: "GPT-5"
  },

  apple: {
    input: "fresh apple",
    description: "🍎 Fruit with natural sugars - fructose, glucose",
    lookupMode: "GPT-5"
  },

  wine: {
    input: "red wine",
    description: "🍷 Alcoholic beverage - ethanol, water, organic acids",
    lookupMode: "GPT-5"
  },

  egg: {
    input: "egg",
    description: "🥚 Protein-rich food - various amino acids, lipids",
    lookupMode: "GPT-5"
  }
};

/**
 * Checks if server is ready to accept requests
 */
async function checkServerHealth() {
  try {
    const response = await fetch(`${SERVER_URL}/health`, {
      timeout: 3000
    });
    return response.ok;
  } catch (error) {
    return false;
  }
}

/**
 * Triggers molecular analysis via HTTP request to live server
 * @param {Object} testCase - Test case configuration
 */
async function injectAnalysis(testCase) {
  console.log(`🔥 Injecting: ${testCase.input}`);
  console.log(`📋 ${testCase.description}`);
  
  try {
    const startTime = Date.now();
    
    // Make request to live server - same as frontend would do
    const response = await fetch(`${SERVER_URL}/structuralize`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify({
        object: testCase.input,
        lookupMode: testCase.lookupMode
      }),
      timeout: 30000 // 30 second timeout
    });
    
    if (!response.ok) {
      throw new Error(`Server returned ${response.status}: ${response.statusText}`);
    }
    
    const result = await response.json();
    const duration = Date.now() - startTime;
    
    console.log(`✅ Analysis completed in ${duration}ms`);
    console.log(`📊 Found ${result.chemicals?.length || 0} chemicals:`);
    
    if (result.chemicals && result.chemicals.length > 0) {
      result.chemicals.slice(0, 5).forEach((chem, i) => {
        const status = chem.status === 'ok' ? '✓' : '⚠️';
        console.log(`   ${i+1}. ${chem.name} ${status}`);
      });
      
      if (result.chemicals.length > 5) {
        console.log(`   ... and ${result.chemicals.length - 5} more molecules`);
      }
    }
    
    console.log(`🎯 Results are now live in the web app!`);
    return result;
    
  } catch (error) {
    console.error(`❌ Injection failed: ${error.message}`);
    
    if (error.message.includes('ECONNREFUSED')) {
      console.error(`💡 Is the server running? Try: ./run start`);
    }
    
    throw error;
  }
}

/**
 * Waits for server to become ready
 */
async function waitForServer(maxWaitMs = 15000) {
  const startTime = Date.now();
  
  while (Date.now() - startTime < maxWaitMs) {
    if (await checkServerHealth()) {
      return true;
    }
    
    console.log('⏳ Waiting for server...');
    await new Promise(resolve => setTimeout(resolve, 1000));
  }
  
  return false;
}

/**
 * Main execution function
 */
async function main() {
  const testName = process.argv[2] || 'coffee';
  
  console.log('🚀 Test Injection Tool');
  console.log('=====================');
  console.log(`Target: ${SERVER_URL}`);
  console.log(`Test case: ${testName}`);
  console.log('');
  
  // Check if test case exists
  if (!TEST_CASES[testName]) {
    console.error(`❌ Unknown test case: ${testName}`);
    console.log(`📋 Available tests: ${Object.keys(TEST_CASES).join(', ')}`);
    process.exit(1);
  }
  
  const testCase = TEST_CASES[testName];
  
  try {
    // Ensure server is ready
    console.log('🔍 Checking server status...');
    const serverReady = await waitForServer();
    
    if (!serverReady) {
      throw new Error('Server is not responding. Make sure it\'s running.');
    }
    
    console.log('✅ Server is ready!');
    console.log('');
    
    // Inject the test
    await injectAnalysis(testCase);
    
    console.log('');
    console.log('🎊 Injection complete!');
    console.log('   → Check your browser for live molecular visualization');
    console.log('   → Results should appear in the web app immediately');
    
  } catch (error) {
    console.error('❌ Test injection failed:', error.message);
    process.exit(1);
  }
}

// Export for use in other scripts
module.exports = {
  TEST_CASES,
  injectAnalysis,
  checkServerHealth,
  waitForServer
};

// Run if called directly
if (require.main === module) {
  main().catch(console.error);
}
