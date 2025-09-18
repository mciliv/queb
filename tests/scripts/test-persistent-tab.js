#!/usr/bin/env node

const ChromeTabManager = require('./test/utils/chrome-tab-manager');

async function runPersistentTabTest() {
  console.log('🎯 Persistent Tab Molecular Testing');
  console.log('===================================');
  console.log('This keeps using the SAME Chrome tab for all tests!\n');

  const tabManager = new ChromeTabManager();
  
  const testCases = [
    { name: 'Water', input: 'water', emoji: '💧' },
    { name: 'Ethanol', input: 'ethanol', emoji: '🧪' },
    { name: 'Coffee', input: 'black coffee', emoji: '☕' },
    { name: 'Wine', input: 'red wine', emoji: '🍷' },
    { name: 'Salt', input: 'sodium chloride', emoji: '🧂' }
  ];

  try {
    // Get or create the persistent tab
    console.log('🔍 Looking for existing molecular testing tab...');
    const page = await tabManager.getOrCreateTab();
    
    const tabInfo = await tabManager.getTabInfo();
    console.log(`📄 Using tab: ${tabInfo.title} (${tabInfo.url})\n`);

    // Run tests in the same persistent tab
    for (let i = 0; i < testCases.length; i++) {
      const test = testCases[i];
      console.log(`${test.emoji} Test ${i+1}/${testCases.length}: ${test.name}`);
      console.log(`   Input: "${test.input}"`);

      // Clear and type input
      console.log(`   ⌨️  Typing "${test.input}"...`);
      await tabManager.typeInput(test.input);

      // Take before screenshot
      const beforePath = await tabManager.takeScreenshot(`${test.name.toLowerCase()}_before`);
      console.log(`   📸 Before: ${beforePath}`);

      // Trigger analysis
      console.log('   🚀 Triggering molecular analysis...');
      await tabManager.triggerAnalysis();

      // Wait for analysis
      console.log('   ⏳ Waiting for molecular visualization...');
      await new Promise(resolve => setTimeout(resolve, 8000));

      // Take after screenshot
      const afterPath = await tabManager.takeScreenshot(`${test.name.toLowerCase()}_result`);
      console.log(`   📸 Result: ${afterPath}`);

      console.log(`   ✅ ${test.name} test complete\n`);

      // Brief pause for observation
      if (i < testCases.length - 1) {
        console.log('   ⏸️  Pausing 3 seconds...');
        await new Promise(resolve => setTimeout(resolve, 3000));
      }
    }

    console.log('🎯 All tests complete!');
    console.log('💡 The same Chrome tab was used for all tests');
    console.log('🔍 You can continue using this tab manually\n');

    // Show final tab info
    const finalTabInfo = await tabManager.getTabInfo();
    console.log(`📄 Final tab state: ${finalTabInfo.title}`);
    console.log(`🌐 URL: ${finalTabInfo.url}`);

    // Cleanup (but keep browser open)
    await tabManager.cleanup();

  } catch (error) {
    console.error('❌ Error during persistent tab testing:', error.message);
    await tabManager.cleanup();
  }
}

// Show usage info
console.log('💡 TIP: If Chrome is already running, this will reuse the existing instance!');
console.log('🔄 Run this script multiple times to keep using the same tab\n');

// Run the test
runPersistentTabTest().catch(console.error);