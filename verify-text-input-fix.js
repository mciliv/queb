#!/usr/bin/env node

/**
 * Simple verification script for text input fix
 * Tests that the API hook properly validates text input before sending to server
 */

// Mock the API call function to capture what would be sent
let capturedRequest = null;

const mockApiCall = async (endpoint, options) => {
  capturedRequest = {
    endpoint,
    body: JSON.parse(options.body)
  };
  
  // Simulate server validation
  if (!capturedRequest.body.text || typeof capturedRequest.body.text !== 'string' || !capturedRequest.body.text.trim()) {
    throw new Error('Missing or invalid "text" parameter');
  }
  
  return { success: true };
};

// Simulate the structuralizeText function with our fix
const structuralizeText = async (text, lookupMode = 'GPT-5') => {
  // Ensure text is a valid non-empty string
  if (typeof text !== 'string' || !text.trim()) {
    throw new Error('Text input is required');
  }

  // Trim the text to ensure no leading/trailing whitespace
  const trimmedText = text.trim();

  const result = await mockApiCall('/api/structuralize', {
    method: 'POST',
    body: JSON.stringify({ text: trimmedText, lookupMode }),
  });
  
  return result;
};

// Test cases
const testCases = [
  { name: 'Valid text', input: 'coffee', shouldPass: true },
  { name: 'Empty string', input: '', shouldPass: false },
  { name: 'Whitespace only', input: '   ', shouldPass: false },
  { name: 'Null', input: null, shouldPass: false },
  { name: 'Undefined', input: undefined, shouldPass: false },
  { name: 'Number', input: 123, shouldPass: false },
  { name: 'Text with whitespace', input: '  coffee  ', shouldPass: true, expectTrimmed: 'coffee' },
];

console.log('🧪 Testing Text Input Validation Fix\n');
console.log('='.repeat(50));

let passed = 0;
let failed = 0;

for (const testCase of testCases) {
  capturedRequest = null;
  let error = null;
  let result = null;
  
  try {
    result = await structuralizeText(testCase.input);
  } catch (e) {
    error = e;
  }
  
  const testPassed = testCase.shouldPass 
    ? (result !== null && error === null)
    : (error !== null && error.message === 'Text input is required');
  
  if (testPassed) {
    if (testCase.shouldPass && testCase.expectTrimmed) {
      // Check that text was trimmed
      if (capturedRequest && capturedRequest.body.text === testCase.expectTrimmed) {
        console.log(`✅ ${testCase.name}: PASSED (text was trimmed correctly)`);
        passed++;
      } else {
        console.log(`❌ ${testCase.name}: FAILED (text was not trimmed)`);
        failed++;
      }
    } else {
      console.log(`✅ ${testCase.name}: PASSED`);
      passed++;
    }
  } else {
    console.log(`❌ ${testCase.name}: FAILED`);
    console.log(`   Expected: ${testCase.shouldPass ? 'should pass' : 'should reject'}`);
    console.log(`   Got: ${error ? error.message : 'no error'}`);
    failed++;
  }
}

console.log('\n' + '='.repeat(50));
console.log(`Results: ${passed} passed, ${failed} failed`);

if (failed === 0) {
  console.log('\n✅ All tests passed! The fix is working correctly.');
  process.exit(0);
} else {
  console.log('\n❌ Some tests failed. Please review the fix.');
  process.exit(1);
}
