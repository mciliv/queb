#!/usr/bin/env node

// Environment validation script to prevent configuration issues
const fs = require('fs');
const path = require('path');

function validateEnvironment() {
  // Find queb project root (directory containing package.json)
  const findProjectRoot = () => {
    let currentDir = __dirname;
    while (currentDir !== path.dirname(currentDir)) {
      if (fs.existsSync(path.join(currentDir, 'package.json'))) {
        return currentDir;
      }
      currentDir = path.dirname(currentDir);
    }
    return path.resolve(__dirname, '..'); // Fallback: go up from scripts to queb
  };
  const projectRoot = findProjectRoot();
  const envPath = path.join(projectRoot, '.env');
  
  if (!fs.existsSync(envPath)) {
    console.error('❌ .env file not found');
    process.exit(1);
  }
  
  const envContent = fs.readFileSync(envPath, 'utf8');
  const lines = envContent.split('\n');
  
  let hasOpenAIModel = false;
  let hasOpenAIKey = false;
  let modelValue = null;
  
  for (const line of lines) {
    const trimmed = line.trim();
    
    // Skip comments and empty lines
    if (trimmed.startsWith('#') || trimmed === '') continue;
    
    if (trimmed.startsWith('OPENAI_MODEL=')) {
      hasOpenAIModel = true;
      modelValue = trimmed.split('=')[1];
    }
    
    if (trimmed.startsWith('OPENAI_API_KEY=')) {
      hasOpenAIKey = true;
    }
  }
  
  // Validation checks
  const errors = [];
  
  if (!hasOpenAIKey) {
    errors.push('Missing OPENAI_API_KEY in .env file');
  }
  
  if (!hasOpenAIModel) {
    errors.push('Missing OPENAI_MODEL in .env file');
  } else if (modelValue && !['gpt-4o', 'gpt-4', 'gpt-3.5-turbo'].includes(modelValue)) {
    errors.push(`Invalid OPENAI_MODEL: ${modelValue}. Use gpt-4o, gpt-4, or gpt-3.5-turbo`);
  }
  
  if (errors.length > 0) {
    console.error('❌ Environment validation failed:');
    errors.forEach(error => console.error(`   - ${error}`));
    console.error('\n💡 Check OpenAI documentation: https://platform.openai.com/docs/models');
    process.exit(1);
  }

  // Silent on success - only report errors
}

if (require.main === module) {
  validateEnvironment();
}

module.exports = { validateEnvironment };






