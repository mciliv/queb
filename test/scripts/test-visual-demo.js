#!/usr/bin/env node

/**
 * Visual Demo Script - Slow, Interactive Molecular Testing
 * Run this to see the molecular visualization tests in action
 */

const puppeteer = require('puppeteer');

// Test cases to demonstrate
const demoTests = [
  { name: '💧 Water', input: 'water', description: 'Simple molecule test' },
  { name: '🧪 Ethanol', input: 'ethanol', description: 'Basic organic compound' },
  { name: '🍷 Red Wine', input: 'red wine', description: 'Complex beverage analysis' },
  { name: '☕ Coffee', input: 'black coffee', description: 'Caffeine detection test' }
];

async function runVisualDemo() {
  console.log('🎬 Starting Visual Molecular Testing Demo');
  console.log('   This will be slow so you can watch each step\n');

  // Launch browser with visible window
  console.log('🚀 Launching browser...');
  const browser = await puppeteer.launch({
    headless: false,
    defaultViewport: null, // Use full window
    args: [
      '--start-maximized',
      '--no-sandbox',
      '--disable-web-security'
    ],
    slowMo: 500 // 500ms delay between actions
  });

  const page = await browser.newPage();
  
  // Log browser console messages
  page.on('console', msg => {
    console.log(`   📄 Browser: ${msg.text()}`);
  });

  try {
    console.log('🌐 Navigating to molecular interface...');
    await page.goto('http://localhost:3001', { waitUntil: 'networkidle0' });
    
    // Wait for user to see the page load
    console.log('⏳ Waiting 3 seconds for you to see the interface...');
    await new Promise(resolve => setTimeout(resolve, 3000));
    
    // Check if app loaded properly
    await page.waitForSelector('#root', { timeout: 10000 });
    console.log('✅ App loaded successfully!\n');

    // Run each test case slowly
    for (let i = 0; i < demoTests.length; i++) {
      const test = demoTests[i];
      console.log(`🧪 Test ${i + 1}/${demoTests.length}: ${test.name}`);
      console.log(`   Description: ${test.description}`);
      console.log(`   Input: "${test.input}"`);

      try {
        // Find the text input
        console.log('   📝 Finding text input...');
        const textInput = await page.waitForSelector('input[type="text"], textarea', { timeout: 5000 });
        
        // Clear and type the test input
        console.log(`   ⌨️  Typing "${test.input}"...`);
        await textInput.click({ clickCount: 3 }); // Select all
        await page.keyboard.type(test.input, { delay: 100 });
        
        // Wait so you can see what was typed
        await new Promise(resolve => setTimeout(resolve, 1500));
        
        // Find and click analyze button
        console.log('   🔍 Looking for analyze button...');
        const buttons = await page.$$('button');
        let analyzeButton = null;
        
        for (const button of buttons) {
          const text = await button.evaluate(el => el.textContent?.toLowerCase() || '');
          if (text.includes('analyze') || text.includes('submit')) {
            analyzeButton = button;
            break;
          }
        }
        
        if (!analyzeButton && buttons.length > 0) {
          analyzeButton = buttons[buttons.length - 1]; // Try last button
        }
        
        if (analyzeButton) {
          console.log('   🚀 Clicking analyze button...');
          await analyzeButton.click();
          
          console.log('   ⏳ Waiting for analysis (you should see the molecular visualization appear)...');
          
          // Wait longer to see results or timeout
          try {
            await page.waitForSelector('.results-section, .molecular-results, .test-results', { timeout: 15000 });
            console.log('   ✅ Results appeared!');
          } catch (e) {
            console.log('   ⚠️  No results section found, but analysis was triggered');
          }
          
          // Wait so you can see the results
          console.log('   👀 Waiting 5 seconds for you to examine the results...');
          await new Promise(resolve => setTimeout(resolve, 5000));
          
        } else {
          console.log('   ❌ Could not find analyze button');
        }
        
      } catch (error) {
        console.log(`   ❌ Error in test: ${error.message}`);
      }
      
      console.log(''); // Empty line between tests
    }

    // Inject fake molecular data as a demo
    console.log('🧪 BONUS: Injecting fake molecular visualization...');
    await page.evaluate(() => {
      const root = document.getElementById('root');
      if (root) {
        const demoDiv = document.createElement('div');
        demoDiv.innerHTML = `
          <div style="position: fixed; top: 20px; right: 20px; background: rgba(0,0,0,0.8); color: white; padding: 20px; border-radius: 10px; z-index: 9999; max-width: 300px;">
            <h3>🧪 Demo Molecular Results</h3>
            <div style="margin: 10px 0;">
              <strong>Water (H₂O)</strong><br>
              <div style="width: 100px; height: 100px; background: #4CAF50; border-radius: 50%; margin: 10px 0; display: flex; align-items: center; justify-content: center; color: white;">
                💧 3D Model
              </div>
            </div>
            <div style="margin: 10px 0;">
              <strong>Ethanol (C₂H₆O)</strong><br>
              <div style="width: 100px; height: 100px; background: #FF9800; border-radius: 50%; margin: 10px 0; display: flex; align-items: center; justify-content: center; color: white;">
                🧪 3D Model
              </div>
            </div>
          </div>
        `;
        root.appendChild(demoDiv);
        console.log('✅ Fake molecular visualization injected!');
      }
    });

    console.log('\n🎉 Demo completed! The browser will stay open for 30 seconds so you can explore...');
    console.log('   You can:');
    console.log('   • Try typing your own objects in the input field');
    console.log('   • Look for the 🧪 test panel button');
    console.log('   • Examine the fake molecular results on the right');
    
    // Keep browser open for exploration
    await new Promise(resolve => setTimeout(resolve, 30000));

  } catch (error) {
    console.error('❌ Demo failed:', error.message);
  } finally {
    console.log('\n👋 Closing browser...');
    await browser.close();
    console.log('✅ Demo finished!');
  }
}

// Run the demo
if (require.main === module) {
  runVisualDemo().catch(error => {
    console.error('❌ Fatal error:', error);
    process.exit(1);
  });
}