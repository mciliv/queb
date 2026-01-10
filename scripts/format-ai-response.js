#!/usr/bin/env node

/**
 * Utility script to format AI responses in a readable way
 * Usage: node scripts/format-ai-response.js '{"object": "[object Object]", "chemicals": []}'
 */

const response = process.argv[2];

if (!response) {
  console.log('Usage: node scripts/format-ai-response.js <json-string>');
  console.log('Example: node scripts/format-ai-response.js \'{"object": "coffee", "chemicals": []}\'');
  process.exit(1);
}

try {
  // First try parsing as-is
  let parsed = JSON.parse(response);
  console.log('📋 Formatted AI Response:');
  console.log('═'.repeat(50));
  console.log(JSON.stringify(parsed, null, 2));
  console.log('═'.repeat(50));
} catch (firstError) {
  try {
    // Try parsing with escaped characters unescaped
    const unescaped = response.replace(/\\n/g, '\n').replace(/\\"/g, '"').replace(/\\'/g, "'");
    const parsed = JSON.parse(unescaped);
    console.log('📋 Formatted AI Response (unescaped):');
    console.log('═'.repeat(50));
    console.log(JSON.stringify(parsed, null, 2));
    console.log('═'.repeat(50));
  } catch (secondError) {
    console.error('❌ Failed to parse JSON with both methods:');
    console.error('Direct parse error:', firstError.message);
    console.error('Unescaped parse error:', secondError.message);
    console.log('Raw input:', response);
  }
}
