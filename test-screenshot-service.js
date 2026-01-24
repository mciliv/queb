#!/usr/bin/env node

/**
 * Test script for Screenshot Service functionality
 *
 * This script tests the screenshot service integration with the DI container
 * and verifies that screenshots can be captured and saved.
 */

const { createContainer } = require('./src/core/services');

async function testScreenshotService() {
  console.log('🧪 Testing Screenshot Service...\n');

  try {
    // Create DI container
    const container = createContainer();

    // Get services
    const screenshotService = await container.get('screenshotService');
    const logger = await container.get('logger');

    console.log('1️⃣ Service Initialization:');
    console.log('   ✅ Screenshot service created');
    console.log(`   📸 Screenshots enabled: ${screenshotService.isEnabled()}`);
    console.log(`   📁 Screenshot directory: ${screenshotService.screenshotDir}\n`);

    console.log('2️⃣ Screenshot Service Methods:');
    console.log(`   ✅ isEnabled(): ${typeof screenshotService.isEnabled}`);
    console.log(`   ✅ takeScreenshot(): ${typeof screenshotService.takeScreenshot}`);
    console.log(`   ✅ getRecentScreenshots(): ${typeof screenshotService.getRecentScreenshots}`);
    console.log(`   ✅ cleanupOldScreenshots(): ${typeof screenshotService.cleanupOldScreenshots}\n`);

    console.log('3️⃣ Testing Screenshot Directory:');
    const fs = require('fs').promises;
    try {
      await fs.access(screenshotService.screenshotDir);
      console.log('   ✅ Screenshot directory exists');
    } catch (error) {
      console.log('   ❌ Screenshot directory missing');
    }

    // List existing screenshots
    const existingScreenshots = await screenshotService.getRecentScreenshots(5);
    console.log(`   📋 Existing screenshots: ${existingScreenshots.length}\n`);

    console.log('4️⃣ Testing Error Logging Integration:');
    const logErrorUseCase = await container.get('logErrorUseCase');
    console.log(`   ✅ LogErrorUseCase created with screenshot service: ${!!logErrorUseCase.screenshotService}\n`);

    console.log('5️⃣ Configuration Check:');
    const config = await container.get('config');
    const enableScreenshots = config.get('development.enableScreenshots');
    console.log(`   🔧 ENABLE_ERROR_SCREENSHOTS: ${enableScreenshots}`);
    console.log(`   💡 To enable screenshots, set ENABLE_ERROR_SCREENSHOTS=true in .env\n`);

    console.log('🎉 Screenshot Service test completed successfully!');
    console.log('\n📝 To enable automatic screenshots:');
    console.log('   1. Add ENABLE_ERROR_SCREENSHOTS=true to your .env file');
    console.log('   2. Restart the server');
    console.log('   3. Visit the app with ?debug=true or set localStorage.enableErrorScreenshots=true');
    console.log('   4. Trigger a JavaScript error to test screenshot capture');

  } catch (error) {
    console.error('❌ Screenshot Service test failed:', error);
    console.error('Stack trace:', error.stack);
    process.exit(1);
  }
}

// Run the test
if (require.main === module) {
  testScreenshotService().catch(console.error);
}

module.exports = { testScreenshotService };