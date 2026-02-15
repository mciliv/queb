#!/usr/bin/env node
/**
 * Visual Test Runner - Automated molecular visualization testing
 * 
 * This script runs chemical analysis on predefined test inputs,
 * generates molecular data, and automatically opens the web app
 * to display the 3D structures for visual verification.
 * 
 * Usage: node scripts/run-visual-test.js [test-name]
 */

const path = require('path');
const fs = require('fs');
const { spawn } = require('child_process');
const appPlacement = require('../config/app-placement');

// Predefined test cases with specific molecular expectations
const TEST_CASES = {
  coffee: {
    input: "coffee",
    description: "Coffee contains caffeine, water, and various organic compounds",
    expectedMolecules: ["caffeine", "water", "chlorogenic acid"],
    lookupMode: "GPT-5"
  },
  
  saltwater: {
    input: "salt water", 
    description: "Simple ionic solution for testing ionic compound visualization",
    expectedMolecules: ["sodium chloride", "water"],
    lookupMode: "GPT-5"
  },
  
  aspirin: {
    input: "aspirin tablet",
    description: "Common medication for testing pharmaceutical compounds", 
    expectedMolecules: ["acetylsalicylic acid", "microcrystalline cellulose"],
    lookupMode: "GPT-5"
  },
  
  apple: {
    input: "fresh apple",
    description: "Fruit with natural sugars and organic acids",
    expectedMolecules: ["fructose", "glucose", "malic acid", "water"],
    lookupMode: "GPT-5"
  },
  
  wine: {
    input: "red wine",
    description: "Alcoholic beverage with complex organic compounds",
    expectedMolecules: ["ethanol", "water", "resveratrol", "tartaric acid"],
    lookupMode: "GPT-5"
  }
};

/**
 * Creates mock analysis results for testing when AI is unavailable
 * @param {Object} testCase - Test case configuration
 * @returns {Object} Mock analysis results
 */
function createMockAnalysis(testCase) {
  const mockData = {
    coffee: {
      object: 'coffee',
      chemicals: [
        { name: 'caffeine', sdfPath: '/sdf_files/caffeine_mock.sdf', status: 'ok' },
        { name: 'water', sdfPath: '/sdf_files/water_mock.sdf', status: 'ok' },
        { name: 'chlorogenic acid', sdfPath: '/sdf_files/chlorogenic_acid_mock.sdf', status: 'ok' }
      ]
    },
    saltwater: {
      object: 'salt water', 
      chemicals: [
        { name: 'sodium chloride', sdfPath: '/sdf_files/sodium_chloride_mock.sdf', status: 'ok' },
        { name: 'water', sdfPath: '/sdf_files/water_mock.sdf', status: 'ok' }
      ]
    },
    aspirin: {
      object: 'aspirin tablet',
      chemicals: [
        { name: 'acetylsalicylic acid', sdfPath: '/sdf_files/acetylsalicylic_acid_mock.sdf', status: 'ok' }
      ]
    },
    apple: {
      object: 'fresh apple',
      chemicals: [
        { name: 'fructose', sdfPath: '/sdf_files/fructose_mock.sdf', status: 'ok' },
        { name: 'glucose', sdfPath: '/sdf_files/glucose_mock.sdf', status: 'ok' },
        { name: 'water', sdfPath: '/sdf_files/water_mock.sdf', status: 'ok' }
      ]
    },
    wine: {
      object: 'red wine',
      chemicals: [
        { name: 'ethanol', sdfPath: '/sdf_files/ethanol_mock.sdf', status: 'ok' },
        { name: 'water', sdfPath: '/sdf_files/water_mock.sdf', status: 'ok' }
      ]
    }
  };
  
  return mockData[testCase.input.toLowerCase().replace(/\s+/g, '')] || mockData.coffee;
}

/**
 * Runs chemical analysis on a test case
 * @param {Object} testCase - Test case configuration
 * @returns {Promise<Object>} Analysis results with molecular data
 */
async function runAnalysis(testCase) {
  console.log(`🧪 Analyzing: ${testCase.input}`);
  console.log(`📋 ${testCase.description}`);
  
  const startTime = Date.now();
  
  // Check if we should use pure mock mode (no API calls)
  const useMockOnly = process.env.MOCK_ONLY === 'true' || process.argv.includes('--mock-only');
  
  if (useMockOnly) {
    console.log(`🤖 Using pure mock data (no API calls)`);
    
    // Simulate processing delay
    await new Promise(resolve => setTimeout(resolve, 500));
    
    const mockResult = createMockAnalysis(testCase);
    const duration = Date.now() - startTime;
    
    console.log(`✅ Mock analysis completed in ${duration}ms`);
    console.log(`📊 Found ${mockResult.chemicals?.length || 0} chemicals:`);
    
    mockResult.chemicals?.forEach((chem, i) => {
      console.log(`   ${i+1}. ${chem.name} (mock data)`);
    });
    
    return mockResult;
  }
  
  // Try real analysis first
  try {
    // Import the analysis service
    const { chemicals } = require('../src/server/services/Structuralizer');
    
    const result = await chemicals({
      object: testCase.input,
      lookupMode: testCase.lookupMode,
      testConfig: {
        // Use faster model for testing
        model: 'gpt-3.5-turbo'
      }
    });
    
    const duration = Date.now() - startTime;
    console.log(`✅ Real analysis completed in ${duration}ms`);
    console.log(`📊 Found ${result.chemicals?.length || 0} chemicals:`);
    
    result.chemicals?.forEach((chem, i) => {
      console.log(`   ${i+1}. ${chem.name} (${chem.status})`);
    });
    
    return result;
    
  } catch (error) {
    console.warn(`⚠️ Real analysis failed: ${error.message}`);
    console.log(`🤖 Falling back to mock data for visualization`);
    
    // Use mock data when real analysis fails
    const mockResult = createMockAnalysis(testCase);
    const duration = Date.now() - startTime;
    
    console.log(`✅ Mock fallback completed in ${duration}ms`);
    console.log(`📊 Found ${mockResult.chemicals?.length || 0} chemicals:`);
    
    mockResult.chemicals?.forEach((chem, i) => {
      console.log(`   ${i+1}. ${chem.name} (mock data)`);
    });
    
    return mockResult;
  }
}

/**
 * Saves test results to a file that the web app can load
 * @param {Object} results - Analysis results
 * @param {string} testName - Name of the test case
 */
async function saveTestResults(results, testName) {
  const testDataDir = path.join(__dirname, '../tests/fixtures/visual-test-data');
  
  // Ensure test data directory exists
  if (!fs.existsSync(testDataDir)) {
    fs.mkdirSync(testDataDir, { recursive: true });
  }
  
  const testDataPath = path.join(testDataDir, `${testName}.json`);
  
  // Prepare data for visualization
  const visualizationData = {
    testName,
    timestamp: new Date().toISOString(),
    object: results.object,
    chemicals: results.chemicals?.map(chem => ({
      name: chem.name,
      sdfPath: chem.sdfPath,
      status: chem.status,
      // Convert file path to URL for browser access
      sdfUrl: chem.sdfPath ? chem.sdfPath.replace(/^.*\/sdf_files\//, '/sdf_files/') : null
    })) || [],
    metadata: {
      totalMolecules: results.chemicals?.length || 0,
      successfulSdfs: results.chemicals?.filter(c => c.status === 'ok').length || 0,
      reason: results.reason
    }
  };
  
  fs.writeFileSync(testDataPath, JSON.stringify(visualizationData, null, 2));
  console.log(`💾 Test data saved to: ${testDataPath}`);
  
  return visualizationData;
}

/**
 * Opens the web app with test data loaded
 * @param {string} testName - Name of the test case to load
 */
async function openVisualization(testName) {
  console.log(`🌐 Opening visualization for: ${testName}`);
  
  // Determine the appropriate URL
  const config = appPlacement.urls.desktop;
  const port = config.defaultPort;
  const protocol = config.protocol;
  
  const testUrl = `${protocol}://localhost:${port}?test=${testName}&autoload=true`;
  
  console.log(`🔗 Opening: ${testUrl}`);
  
  // Use system browser if available, otherwise try Chrome
  const isWindows = process.platform === 'win32';
  const isMac = process.platform === 'darwin';
  const isLinux = process.platform === 'linux';
  
  let browserCmd;
  let browserArgs = [testUrl];
  
  if (isMac) {
    browserCmd = 'open';
  } else if (isWindows) {
    browserCmd = 'start';
    browserArgs = ['', testUrl]; // Windows start command needs empty string
  } else if (isLinux) {
    browserCmd = 'xdg-open';
  } else {
    throw new Error('Unsupported platform for browser automation');
  }
  
  return new Promise((resolve, reject) => {
    const browserProcess = spawn(browserCmd, browserArgs, { 
      stdio: 'ignore',
      detached: true
    });
    
    browserProcess.on('error', (error) => {
      console.warn(`⚠️ Browser open failed: ${error.message}`);
      console.log(`📌 Manual URL: ${testUrl}`);
      resolve(); // Don't fail the test, just provide manual URL
    });
    
    browserProcess.on('spawn', () => {
      console.log(`✅ Browser opened successfully`);
      resolve();
    });
    
    // Detach the process so it doesn't block the script
    browserProcess.unref();
  });
}

/**
 * Main execution function
 */
async function main() {
  const testName = process.argv[2] || 'coffee';
  
  if (!TEST_CASES[testName]) {
    console.error(`❌ Unknown test case: ${testName}`);
    console.log(`📋 Available tests: ${Object.keys(TEST_CASES).join(', ')}`);
    process.exit(1);
  }
  
  const testCase = TEST_CASES[testName];
  
  console.log(`🚀 Running visual test: ${testName}`);
  console.log(`🔬 Input: "${testCase.input}"`);
  console.log(`⚙️ Mode: ${testCase.lookupMode}`);
  console.log('');
  
  try {
    // Step 1: Run chemical analysis
    const results = await runAnalysis(testCase);
    
    // Step 2: Save results for web app
    const visualizationData = await saveTestResults(results, testName);
    
    console.log('');
    console.log(`📊 Test Results Summary:`);
    console.log(`   Object: ${visualizationData.object}`);
    console.log(`   Total molecules: ${visualizationData.metadata.totalMolecules}`);
    console.log(`   With 3D structures: ${visualizationData.metadata.successfulSdfs}`);
    console.log(`   Success rate: ${Math.round((visualizationData.metadata.successfulSdfs / visualizationData.metadata.totalMolecules) * 100)}%`);
    console.log('');
    
    // Step 3: Open browser for visualization
    await openVisualization(testName);
    
    console.log('🎯 Visual test setup complete!');
    console.log('   The web app should open with your test data loaded.');
    console.log('   Verify the 3D molecular structures appear correctly.');
    
  } catch (error) {
    console.error(`❌ Visual test failed: ${error.message}`);
    process.exit(1);
  }
}

// Run if called directly
if (require.main === module) {
  main().catch(console.error);
}

module.exports = {
  TEST_CASES,
  runAnalysis,
  saveTestResults,
  openVisualization,
  main
};
