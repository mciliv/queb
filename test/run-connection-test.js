#!/usr/bin/env node

// run-connection-test.js - Run backwards connection test with actual server
const BackwardsConnectionTester = require('./backwards-connection-test');

async function runConnectionTest() {
  console.log('🚀 Starting Connection Test with actual server...');
  
  try {
    // Import the actual server
    const server = require('../backend/api/server');
    
    // Wait for server to be ready
    await new Promise(resolve => setTimeout(resolve, 2000));
    
    const tester = new BackwardsConnectionTester(server);
    const report = await tester.runBackwardsConnectionTest();
    
    console.log('\n🎯 Test completed with health:', report.connectionHealth + '%');
    
    if (report.connectionHealth >= 60) {
      console.log('✅ Connection test passed');
      process.exit(0);
    } else {
      console.log('❌ Connection test failed');
      process.exit(1);
    }
    
  } catch (error) {
    console.error('❌ Failed to run connection test:', error.message);
    process.exit(1);
  }
}

runConnectionTest(); 