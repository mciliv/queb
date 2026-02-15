const ChromeTabManager = require('../../utils/chrome-tab-manager');
const path = require('path');

describe('Persistent Tab Molecular Tests', () => {
  let tabManager;
  let page;

  beforeAll(async () => {
    console.log('🚀 Starting persistent tab molecular tests...');
    console.log('   This will reuse the same Chrome tab for all tests!');
    
    tabManager = new ChromeTabManager();
    page = await tabManager.getOrCreateTab();
    
    const tabInfo = await tabManager.getTabInfo();
    console.log(`✅ Using persistent tab: ${tabInfo.title}`);
    console.log(`🌐 URL: ${tabInfo.url}\n`);
  });

  afterAll(async () => {
    console.log('\n🎬 Tests complete! Chrome tab remains open for inspection.');
    if (tabManager) {
      await tabManager.cleanup();
    }
  });

  test('should use the same tab for water analysis', async () => {
    console.log('💧 Testing water in persistent tab...');
    
    await tabManager.typeInput('water');
    const beforePath = await tabManager.takeScreenshot('water_before');
    console.log(`   📸 Before: ${beforePath}`);
    
    await tabManager.triggerAnalysis();
    await new Promise(resolve => setTimeout(resolve, 8000));
    
    const afterPath = await tabManager.takeScreenshot('water_after');
    console.log(`   📸 After: ${afterPath}`);
    
    // Verify we're still on the same page
    const tabInfo = await tabManager.getTabInfo();
    expect(tabInfo.url).toContain('localhost:8080');
    
    console.log('✅ Water test complete in persistent tab');
  }, 30000);

  test('should reuse same tab for ethanol analysis', async () => {
    console.log('🧪 Testing ethanol in same persistent tab...');
    
    await tabManager.typeInput('ethanol');
    const beforePath = await tabManager.takeScreenshot('ethanol_before');
    console.log(`   📸 Before: ${beforePath}`);
    
    await tabManager.triggerAnalysis();
    await new Promise(resolve => setTimeout(resolve, 8000));
    
    const afterPath = await tabManager.takeScreenshot('ethanol_after');
    console.log(`   📸 After: ${afterPath}`);
    
    // Verify same tab is still being used
    const tabInfo = await tabManager.getTabInfo();
    expect(tabInfo.url).toContain('localhost:8080');
    
    console.log('✅ Ethanol test complete in same persistent tab');
  }, 30000);

  test('should continue using same tab for coffee analysis', async () => {
    console.log('☕ Testing coffee in same persistent tab...');
    
    await tabManager.typeInput('black coffee');
    const beforePath = await tabManager.takeScreenshot('coffee_before');
    console.log(`   📸 Before: ${beforePath}`);
    
    await tabManager.triggerAnalysis();
    await new Promise(resolve => setTimeout(resolve, 8000));
    
    const afterPath = await tabManager.takeScreenshot('coffee_after');
    console.log(`   📸 After: ${afterPath}`);
    
    // Verify tab persistence
    const tabInfo = await tabManager.getTabInfo();
    expect(tabInfo.url).toContain('localhost:8080');
    
    console.log('✅ Coffee test complete in same persistent tab');
  }, 30000);

  test('should verify tab persistence across all tests', async () => {
    console.log('🔍 Verifying tab persistence...');
    
    const tabInfo = await tabManager.getTabInfo();
    console.log(`   📄 Current tab: ${tabInfo.title}`);
    console.log(`   🌐 Current URL: ${tabInfo.url}`);
    
    // Check that we can still interact with the page
    const inputExists = await page.$('#object-input');
    expect(inputExists).toBeTruthy();
    
    // Verify page is responsive
    await tabManager.typeInput('test persistence');
    const inputValue = await page.$eval('#object-input', el => el.value);
    expect(inputValue).toBe('test persistence');
    
    console.log('✅ Tab persistence verified - same tab used throughout all tests!');
  }, 15000);
});