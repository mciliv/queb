const { chemicals, Structuralizer } = require('../src/server/services/Structuralizer');
const testConfigs = require('./test-config-example');

// Example usage of the new test configuration interface
async function testStructuralizer() {
  console.log('Testing Structuralizer with different configurations...\n');
  
  // Test 1: Custom prompt (uses original requestParams.messages)
  console.log('1. Testing with custom prompt (uses original messages):');
  const structuralizer1 = new Structuralizer(null, testConfigs.customPrompt);
  console.log(`   Model: ${structuralizer1.testConfig.model}`);
  console.log(`   Prompt: ${structuralizer1.testConfig.prompt}`);
  console.log(`   Is test mode: ${structuralizer1.isTestMode}\n`);
  
  // Test 2: Hardcoded prompt
  console.log('2. Testing with hardcoded prompt:');
  const structuralizer2 = new Structuralizer(null, testConfigs.hardcodedPrompt);
  console.log(`   Model: ${structuralizer2.testConfig.model}`);
  console.log(`   Prompt: ${structuralizer2.testConfig.prompt}\n`);
  
  // Test 3: Structuralize test
  console.log('3. Testing with structuralize prompt:');
  const structuralizer3 = new Structuralizer(null, testConfigs.structuralizeTest);
  console.log(`   Model: ${structuralizer3.testConfig.model}`);
  console.log(`   Prompt: ${structuralizer3.testConfig.prompt}\n`);
  const example = await chemicals({ object: 'ethanol' });
  console.log('   Example chemicals():', example);
  
  // Test 4: Object detection test
  console.log('4. Testing with object detection prompt:');
  const structuralizer4 = new Structuralizer(null, testConfigs.objectDetectionTest);
  console.log(`   Model: ${structuralizer4.testConfig.model}`);
  console.log(`   Prompt: ${structuralizer4.testConfig.prompt}\n`);
  
  // Test 5: No config (uses original behavior)
  console.log('5. Testing with no config (original behavior):');
  const structuralizer5 = new Structuralizer();
  console.log(`   Model candidates: ${structuralizer5.modelCandidates.join(', ')}`);
  console.log(`   Test config: ${JSON.stringify(structuralizer5.testConfig)}`);
}

// Run the test
if (require.main === module) {
  testStructuralizer().catch(console.error);
}

module.exports = { testStructuralizer };
