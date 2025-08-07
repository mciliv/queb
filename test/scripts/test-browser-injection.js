/**
 * Browser injection script for testing molecular visualization
 * 
 * Instructions:
 * 1. Open http://localhost:3001 
 * 2. Open browser console (F12)
 * 3. Copy and paste this entire script
 * 4. Watch the molecular visualization appear!
 */

// Test data injection function
async function injectMolecularTest(testName = "Water Test") {
  console.log(`🧪 Running ${testName}...`);
  
  // Find the text input and submit button
  const textInput = document.querySelector('input[type="text"]') || 
                   document.querySelector('textarea') ||
                   document.querySelector('#object-input') ||
                   document.querySelector('[placeholder*="object"]');
  
  const submitButton = document.querySelector('button[type="submit"]') ||
                      document.querySelector('button:contains("Analyze")') ||
                      document.querySelector('.analyze-button') ||
                      document.querySelector('button');

  if (!textInput) {
    console.error('❌ Could not find text input field');
    return;
  }

  if (!submitButton) {
    console.error('❌ Could not find submit button');
    return;
  }

  // Test cases
  const testCases = {
    "Water Test": "water",
    "Wine Test": "red wine", 
    "Coffee Test": "black coffee",
    "Ethanol Test": "ethanol"
  };

  const testInput = testCases[testName] || "water";
  
  console.log(`📝 Injecting "${testInput}" into input field...`);
  
  // Set the input value
  textInput.value = testInput;
  textInput.dispatchEvent(new Event('input', { bubbles: true }));
  textInput.dispatchEvent(new Event('change', { bubbles: true }));
  
  console.log(`🚀 Triggering analysis...`);
  
  // Click the submit button
  submitButton.click();
  
  console.log(`⏳ Analysis started - watch for molecular visualization results!`);
}

// Quick test function for all cases
async function runAllTests() {
  const tests = ["Water Test", "Wine Test", "Coffee Test", "Ethanol Test"];
  
  for (const test of tests) {
    console.log(`\n🔬 Running ${test}...`);
    await injectMolecularTest(test);
    
    // Wait 6 seconds between tests to see results
    await new Promise(resolve => setTimeout(resolve, 6000));
  }
  
  console.log('\n✅ All tests completed!');
}

// Make functions available globally
window.injectMolecularTest = injectMolecularTest;
window.runAllTests = runAllTests;

console.log(`
🧪 MOLECULAR TEST INJECTION READY!

Available commands:
1. injectMolecularTest("Water Test")    - Test single case
2. injectMolecularTest("Wine Test")     - Test complex case  
3. runAllTests()                        - Run all test cases

Just type any of these commands in this console!
`);

// Auto-run a quick test
console.log('🚀 Auto-running Water Test in 2 seconds...');
setTimeout(() => injectMolecularTest("Water Test"), 2000);