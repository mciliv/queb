#!/usr/bin/env node

const VersionManager = require('../utils/version-manager');
const { getAllVersions } = require('../support/multi-version-config');

async function runParallelVersionTests() {
  const manager = new VersionManager();
  
  try {
    console.log('🚀 Starting parallel version testing...');
    
    // Test current and previous versions simultaneously
    const versions = ['main', 'previous'];
    
    console.log(`📋 Testing versions: ${versions.join(', ')}`);
    
    const results = await manager.runParallelTests(versions, 'molecular-basic');
    
    console.log('\n📊 Test Results:');
    for (const [version, result] of results) {
      if (result.error) {
        console.log(`❌ ${version}: ${result.error}`);
      } else {
        console.log(`✅ ${version}: Found ${result.result.moleculesFound} molecules`);
      }
    }
    
    console.log('\n🔄 Versions still running for manual inspection');
    console.log('💡 Use Ctrl+C to cleanup all versions');
    
    // Keep versions running for inspection
    process.on('SIGINT', async () => {
      console.log('\n🧹 Cleaning up all versions...');
      await manager.cleanup();
      process.exit(0);
    });
    
    // Keep process alive
    await new Promise(() => {});
    
  } catch (error) {
    console.error('❌ Parallel testing failed:', error.message);
    await manager.cleanup();
    process.exit(1);
  }
}

// Allow running specific versions
const args = process.argv.slice(2);
if (args.length > 0) {
  const customVersions = args;
  console.log(`🎯 Custom versions specified: ${customVersions.join(', ')}`);
}

runParallelVersionTests();
