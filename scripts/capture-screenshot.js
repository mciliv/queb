#!/usr/bin/env node

/**
 * CLI utility for capturing screenshots for LLM analysis
 * Usage:
 *   node scripts/capture-screenshot.js                    # Basic app screenshot
 *   node scripts/capture-screenshot.js --input "water"    # Screenshot with input
 *   node scripts/capture-screenshot.js --analysis "water" # Screenshot after analysis
 *   node scripts/capture-screenshot.js --list             # List all screenshots
 */

const ScreenshotService = require('../backend/services/screenshot-service');

async function main() {
  const args = process.argv.slice(2);
  const screenshotService = new ScreenshotService();

  try {
    if (args.includes('--list')) {
      // List all screenshots
      console.log('📋 Listing all screenshots...');
      const screenshots = await screenshotService.listScreenshots();
      
      if (screenshots.length === 0) {
        console.log('📭 No screenshots found');
        return;
      }

      console.log(`\n📸 Found ${screenshots.length} screenshots:`);
      screenshots.forEach((screenshot, index) => {
        console.log(`${index + 1}. ${screenshot.filename}`);
        console.log(`   🔗 http://localhost:3000${screenshot.url}`);
      });

    } else if (args.includes('--input')) {
      // Screenshot with input
      const inputIndex = args.indexOf('--input');
      const inputText = args[inputIndex + 1];
      
      if (!inputText) {
        console.error('❌ Please provide input text: --input "water"');
        process.exit(1);
      }

      console.log(`📝 Capturing screenshot with input: "${inputText}"`);
      const result = await screenshotService.captureWithInput(inputText);
      
      console.log('✅ Screenshot captured successfully!');
      console.log(`📁 File: ${result.filename}`);
      console.log(`🔗 URL: http://localhost:3000${result.url}`);

    } else if (args.includes('--analysis')) {
      // Screenshot after analysis
      const analysisIndex = args.indexOf('--analysis');
      const inputText = args[analysisIndex + 1];
      
      if (!inputText) {
        console.error('❌ Please provide input text: --analysis "water"');
        process.exit(1);
      }

      console.log(`🔬 Capturing analysis screenshot for: "${inputText}"`);
      const result = await screenshotService.captureAnalysis(inputText);
      
      console.log('✅ Analysis screenshot captured successfully!');
      console.log(`📁 File: ${result.filename}`);
      console.log(`🔗 URL: http://localhost:3000${result.url}`);

    } else {
      // Basic app screenshot
      console.log('📸 Capturing basic app screenshot...');
      const result = await screenshotService.captureApp();
      
      console.log('✅ Screenshot captured successfully!');
      console.log(`📁 File: ${result.filename}`);
      console.log(`🔗 URL: http://localhost:3000${result.url}`);
    }

  } catch (error) {
    console.error('❌ Screenshot capture failed:', error.message);
    
    if (error.message.includes('ECONNREFUSED')) {
      console.log('💡 Make sure the frontend is running on http://localhost:3000');
    }
    
    process.exit(1);
  } finally {
    await screenshotService.close();
  }
}

// Show usage if --help
if (process.argv.includes('--help') || process.argv.includes('-h')) {
  console.log(`
📸 Screenshot Capture Tool for LLM Analysis

Usage:
  node scripts/capture-screenshot.js                    # Basic app screenshot
  node scripts/capture-screenshot.js --input "water"    # Screenshot with input
  node scripts/capture-screenshot.js --analysis "water" # Screenshot after analysis
  node scripts/capture-screenshot.js --list             # List all screenshots
  node scripts/capture-screenshot.js --help             # Show this help

Examples:
  node scripts/capture-screenshot.js --input "caffeine"
  node scripts/capture-screenshot.js --analysis "aspirin"
  node scripts/capture-screenshot.js --list

Screenshots are saved to: /screenshots/
LLM-accessible URLs: http://localhost:3000/api/screenshot/[filename]
`);
  process.exit(0);
}

main();
