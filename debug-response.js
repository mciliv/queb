#!/usr/bin/env node

/**
 * Debug utility for formatting AI responses viewed in Cursor debugger
 * Copy-paste the variable value from Cursor's debugger and run this script
 */

// Replace this with your variable content from Cursor debugger
// Copy-paste the entire object/array/string from Cursor's Variables panel

const responseFromCursor = {
  // PASTE YOUR VARIABLE HERE
  object: "[object Object]",
  chemicals: [],
  reason: "The provided label appears to be a generic JavaScript placeholder rather than a specific item; without a concrete object name or description, the chemical composition cannot be determined."
};

// Alternative: If you have a JSON string variable
// const jsonString = '{"object": "[object Object]", "chemicals": [], "reason": "..."}';
// const responseFromCursor = JSON.parse(jsonString);

// Utility functions for different variable types
function formatVariable(variable, label = 'Variable') {
  console.log(`\n📋 ${label}:`);
  console.log('─'.repeat(40));

  if (typeof variable === 'string') {
    // Check if it's a JSON string
    try {
      const parsed = JSON.parse(variable);
      console.log('String contains JSON - parsed version:');
      console.log(JSON.stringify(parsed, null, 2));
      console.log('\nRaw string:');
      console.log(variable);
    } catch {
      console.log('Raw string:');
      console.log(variable);
    }
  } else if (typeof variable === 'object') {
    console.log('Object/Array:');
    console.log(variable);
    console.log('\nFormatted JSON:');
    console.log(JSON.stringify(variable, null, 2));
  } else {
    console.log(`${typeof variable}:`, variable);
  }
}

// Format the main variable
formatVariable(responseFromCursor, 'AI Response from Cursor');

// Example for OpenAI API response structure
if (responseFromCursor.choices?.[0]?.message?.content) {
  console.log('\n🔍 Extracted from OpenAI Response:');
  try {
    const content = responseFromCursor.choices[0].message.content;
    const parsed = JSON.parse(content);
    formatVariable(parsed, 'Parsed AI Content');
  } catch (e) {
    console.log('Failed to parse content:', e.message);
    formatVariable(responseFromCursor.choices[0].message.content, 'Raw AI Content');
  }
}
