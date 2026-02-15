#!/usr/bin/env node
/**
 * Test with Server - Complete visual testing workflow
 * 
 * This script ensures the server is running, then executes visual tests
 * with automatic browser opening for result visualization.
 * 
 * Usage: node scripts/test-with-server.js [test-name]
 */

const { spawn } = require('child_process');
const path = require('path');
const fetch = require('node-fetch').default || require('node-fetch');

const SERVER_PORT = 8080;
const SERVER_URL = `http://localhost:${SERVER_PORT}`;
const MAX_SERVER_WAIT = 30000; // 30 seconds max wait for server

/**
 * Checks if the server is running and ready
 * @returns {Promise<boolean>} True if server is accessible
 */
async function isServerRunning() {
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
 * Waits for server to become available
 * @returns {Promise<boolean>} True if server started successfully
 */
async function waitForServer() {
  const startTime = Date.now();
  
  while (Date.now() - startTime < MAX_SERVER_WAIT) {
    if (await isServerRunning()) {
      console.log('✅ Server is ready');
      return true;
    }
    
    console.log('⏳ Waiting for server...');
    await new Promise(resolve => setTimeout(resolve, 2000));
  }
  
  return false;
}

/**
 * Starts the server in background
 * @returns {Promise<ChildProcess>} Server process
 */
async function startServer() {
  console.log('🚀 Starting server...');
  
  const serverProcess = spawn('node', ['src/server/api/server.js'], {
    cwd: path.join(__dirname, '..'),
    stdio: ['ignore', 'inherit', 'inherit'],
    detached: false,
    env: {
      ...process.env
    }
  });
  
  serverProcess.on('error', (error) => {
    console.error('❌ Server start failed:', error.message);
  });
  
  // Give server time to initialize
  await new Promise(resolve => setTimeout(resolve, 3000));
  
  return serverProcess;
}

/**
 * Runs the visual test
 * @param {string} testName - Name of test to run
 * @returns {Promise<void>}
 */
async function runVisualTest(testName) {
  console.log(`🧪 Running visual test: ${testName}`);
  
  return new Promise((resolve, reject) => {
    const testProcess = spawn('node', ['scripts/run-visual-test.js', testName], {
      cwd: path.join(__dirname, '..'),
      stdio: 'inherit'
    });
    
    testProcess.on('close', (code) => {
      if (code === 0) {
        console.log('✅ Visual test completed successfully');
        resolve();
      } else {
        reject(new Error(`Visual test failed with exit code ${code}`));
      }
    });
    
    testProcess.on('error', (error) => {
      reject(new Error(`Visual test process error: ${error.message}`));
    });
  });
}

/**
 * Main execution function
 */
async function main() {
  const testName = process.argv[2] || 'coffee';
  
  console.log('🧪 Visual Testing - One Command Solution');
  console.log('========================================');
  console.log(`Running test: ${testName}`);
  console.log('This will: analyze → generate 3D → open browser');
  console.log('');
  
  let serverProcess = null;
  let serverWasRunning = false;
  
  try {
    // Step 1: Ensure server is running
    console.log('🔍 Checking server status...');
    serverWasRunning = await isServerRunning();
    
    if (!serverWasRunning) {
      console.log('🚀 Starting server (this may take a moment)...');
      serverProcess = await startServer();
      
      const serverReady = await waitForServer();
      if (!serverReady) {
        throw new Error('Server failed to start within timeout period');
      }
    } else {
      console.log('✅ Server already running - proceeding with test');
    }
    
    console.log('');
    
    // Step 2: Run the visual test
    console.log('🧬 Running molecular analysis and visualization...');
    await runVisualTest(testName);
    
    console.log('');
    console.log('🎯 SUCCESS! Visual test complete.');
    console.log('   ✓ Browser should now be open with 3D molecular structures');
    console.log('   ✓ Verify the molecules look correct');
    console.log('   ✓ Test data saved for future reference');
    
    if (serverProcess && !serverWasRunning) {
      console.log('');
      console.log('💡 Server will continue running for additional tests');
      console.log('   Use Ctrl+C to stop the server when done');
    }
    
  } catch (error) {
    console.error('❌ Visual test failed:', error.message);
    
    if (serverProcess && !serverWasRunning) {
      console.log('🛑 Cleaning up server...');
      serverProcess.kill();
    }
    
    process.exit(1);
  }
}

// Run if called directly
if (require.main === module) {
  main().catch(console.error);
}

module.exports = {
  isServerRunning,
  waitForServer,
  startServer,
  runVisualTest,
  main
};
