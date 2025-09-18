#!/usr/bin/env node

// Secure Environment Setup Script
// Helps you choose and configure the right environment library

const fs = require('fs');
const path = require('path');
const { execSync } = require('child_process');

console.log('🔐 Secure Environment Setup');
console.log('===========================\n');

// Check what libraries are available
const libraries = {
  convict: checkLibrary('convict'),
  dotenvx: checkLibrary('@dotenvx/dotenvx'),
  vault: checkLibrary('dotenv-vault'),
  config: checkLibrary('config')
};

function checkLibrary(name) {
  try {
    require.resolve(name);
    return '✅ installed';
  } catch {
    return '❌ not installed';
  }
}

console.log('📦 Available Libraries:');
Object.entries(libraries).forEach(([name, status]) => {
  console.log(`   ${name}: ${status}`);
});
console.log('');

function runDemo(demoName, description, command) {
  console.log(`🧪 ${demoName}`);
  console.log(`   ${description}`);
  console.log(`   Command: ${command}`);

  try {
    const result = execSync(command, { encoding: 'utf8', timeout: 10000 });
    console.log('   ✅ Success!');
    if (result.trim()) console.log(`   Output: ${result.trim()}`);
  } catch (error) {
    console.log('   ❌ Failed:', error.message.split('\n')[0]);
  }
  console.log('');
}

console.log('🚀 Running Library Demos...\n');

// Demo 1: Basic convict validation
if (libraries.convict === '✅ installed') {
  runDemo(
    'Convict Validation',
    'Type-safe configuration with validation',
    'node -e "try { const config = require(\'./config/secure-config\'); console.log(\'Valid config loaded!\'); console.log(\'Port:\', config.get(\'server.port\')); } catch(e) { console.error(\'Error:\', e.message); }"'
  );
}

// Demo 2: DotenvX encryption
if (libraries.dotenvx === '✅ installed') {
  runDemo(
    'DotenvX Encryption',
    'File-level encryption for sensitive data',
    'npx dotenvx --help | head -5'
  );
}

// Demo 3: Vault status
if (libraries.vault === '✅ installed') {
  runDemo(
    'Vault Status',
    'Enterprise secrets management',
    'npx dotenv-vault --help | head -5'
  );
}

console.log('📋 Setup Recommendations:');
console.log('');

const recommendations = [
  {
    scenario: 'Quick & Simple',
    library: 'Basic dotenv',
    setup: 'Already working - just add your API keys to ~/.env'
  },
  {
    scenario: 'Team Development',
    library: 'dotenvx',
    setup: 'npx dotenvx encrypt .env.local && npx dotenvx run -- node your-app.js'
  },
  {
    scenario: 'Enterprise Production',
    library: 'dotenv-vault',
    setup: 'npx dotenv-vault login && npx dotenv-vault push'
  },
  {
    scenario: 'Type Safety + Validation',
    library: 'convict + joi',
    setup: 'node -r ./config/secure-config.js your-app.js'
  },
  {
    scenario: 'Complex Configuration',
    library: 'config (npm)',
    setup: 'Create config/*.json files and require(\'config\')'
  }
];

recommendations.forEach(rec => {
  console.log(`• ${rec.scenario}: ${rec.library}`);
  console.log(`  → ${rec.setup}`);
});

console.log('');
console.log('🔧 Quick Setup Commands:');
console.log('');
console.log('1. For immediate use (basic):');
console.log('   cp .env.local.example .env.local');
console.log('   # Edit .env.local with your actual secrets');
console.log('');
console.log('2. For team collaboration (dotenvx):');
console.log('   npx dotenvx encrypt .env.local');
console.log('   npx dotenvx run -- node your-app.js');
console.log('');
console.log('3. For enterprise (vault):');
console.log('   npx dotenv-vault login');
console.log('   npx dotenv-vault push');
console.log('');
console.log('4. For validation (convict):');
console.log('   node -e "const config = require(\'./config/secure-config\'); console.log(config.get());"');
console.log('');
console.log('📚 Documentation:');
console.log('   • Convict: https://github.com/mozilla/node-convict');
console.log('   • DotenvX: https://dotenvx.com');
console.log('   • Vault: https://www.dotenv.org/docs/vault');
console.log('   • Config: https://github.com/node-config/node-config');

console.log('');
console.log('🎯 Choose the approach that fits your needs!');
console.log('   Start simple, upgrade as your project grows. 🔐');
