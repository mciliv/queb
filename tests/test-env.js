#!/usr/bin/env node

// Test script for the simplified environment loading system
const env = require('../config/env');

console.log('🧪 Testing Simplified Environment Loading\n');

// Test validation
console.log('1. Validation test:');
try {
  env.validate(['NODE_ENV']);
  console.log('✅ All required variables present');
} catch (error) {
  console.log('❌', error.message);
}

// Test specific variables
console.log('\n4. Key environment variables:');
const keyVars = ['NODE_ENV', 'OPENAI_API_KEY', 'PORT', 'DB_HOST'];
keyVars.forEach(key => {
  const value = process.env[key];
  const masked = key.includes('KEY') || key.includes('PASSWORD') || key.includes('SECRET')
    ? (value ? '***' + value.slice(-4) : 'not set')
    : (value || 'not set');
  console.log(`   ${key}: ${masked}`);
});

console.log('\n✅ Environment test complete!');
