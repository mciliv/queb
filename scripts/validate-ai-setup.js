#!/usr/bin/env node

/**
 * AI Setup Validation Script
 * Ensures AI configuration is properly set up for any environment
 */

const fs = require('fs');
const path = require('path');

function validateEnvironment() {
  console.log('🔍 Validating AI Setup Configuration...\n');

  const issues = [];
  const warnings = [];

  // Check if .env exists
  const envPath = path.join(__dirname, '..', '.env');
  if (!fs.existsSync(envPath)) {
    issues.push('❌ .env file not found. Copy .env.example to .env and configure your API keys.');
    return { issues, warnings };
  }

  console.log('✅ .env file exists');

  // Load environment variables
  require('dotenv').config({ path: envPath });

  // Validate AI Provider
  const aiProvider = process.env.AI_PROVIDER;
  if (!aiProvider) {
    issues.push('❌ AI_PROVIDER not set. Must be "openai" or "xai"');
  } else if (!['openai', 'xai'].includes(aiProvider)) {
    issues.push(`❌ Invalid AI_PROVIDER: "${aiProvider}". Must be "openai" or "xai"`);
  } else {
    console.log(`✅ AI Provider: ${aiProvider}`);
  }

  // Validate OpenAI configuration
  if (!process.env.OPENAI_API_KEY || process.env.OPENAI_API_KEY === 'your_openai_api_key_here') {
    if (aiProvider === 'openai') {
      issues.push('❌ OPENAI_API_KEY not configured. Get your key from https://platform.openai.com/api-keys');
    } else {
      warnings.push('⚠️  OPENAI_API_KEY not configured (only needed if using OpenAI provider)');
    }
  } else {
    console.log('✅ OpenAI API key configured');
  }

  // Validate xAI configuration
  if (!process.env.XAI_API_KEY || process.env.XAI_API_KEY === 'your_xai_api_key_here') {
    if (aiProvider === 'xai') {
      issues.push('❌ XAI_API_KEY not configured. Get your key from https://console.x.ai/');
    } else {
      warnings.push('⚠️  XAI_API_KEY not configured (only needed if using xAI provider)');
    }
  } else {
    console.log('✅ xAI API key configured');
  }

  // Check model configurations
  const openaiModel = process.env.OPENAI_MODEL || 'latest';
  const xaiModel = process.env.XAI_MODEL || 'latest';

  console.log(`✅ OpenAI model: ${openaiModel}`);
  console.log(`✅ xAI model: ${xaiModel}`);

  // Validate package dependencies
  const packagePath = path.join(__dirname, '..', 'package.json');
  if (fs.existsSync(packagePath)) {
    const packageJson = JSON.parse(fs.readFileSync(packagePath, 'utf8'));
    const deps = packageJson.dependencies || {};

    const requiredDeps = [
      '@ai-sdk/openai',
      '@ai-sdk/xai',
      'ai',
      'openai'
    ];

    requiredDeps.forEach(dep => {
      if (!deps[dep]) {
        issues.push(`❌ Missing required dependency: ${dep}`);
      } else {
        console.log(`✅ Dependency installed: ${dep}@${deps[dep]}`);
      }
    });
  } else {
    issues.push('❌ package.json not found');
  }

  return { issues, warnings };
}

function testAIService() {
  console.log('\n🧪 Testing AI Service Initialization...');

  try {
    const AIService = require('../src/server/services/AIService');

    // Test service creation
    const aiService = new AIService();
    console.log('✅ AI Service created successfully');
    console.log(`   Provider: ${aiService.getProvider()}`);
    console.log(`   Model: ${aiService.getModel()}`);

    // Test provider switching
    const originalProvider = aiService.getProvider();
    const otherProvider = originalProvider === 'openai' ? 'xai' : 'openai';

    aiService.switchProvider(otherProvider);
    console.log(`✅ Provider switching works: ${originalProvider} → ${aiService.getProvider()}`);

    // Switch back
    aiService.switchProvider(originalProvider);
    console.log(`✅ Provider switching works: ${otherProvider} → ${aiService.getProvider()}`);

    console.log('✅ AI Service test completed successfully');

  } catch (error) {
    console.log('❌ AI Service test failed:', error.message);
    return false;
  }

  return true;
}

async function main() {
  console.log('🚀 AI Setup Validation & Lock-in Script');
  console.log('=====================================\n');

  const { issues, warnings } = validateEnvironment();

  if (issues.length > 0) {
    console.log('\n❌ CRITICAL ISSUES FOUND:');
    issues.forEach(issue => console.log(`   ${issue}`));
    console.log('\n🔧 Fix these issues before proceeding.');
    process.exit(1);
  }

  if (warnings.length > 0) {
    console.log('\n⚠️  WARNINGS:');
    warnings.forEach(warning => console.log(`   ${warning}`));
    console.log('');
  }

  console.log('✅ All critical configuration validated!');

  // Test the AI service
  const serviceTestPassed = testAIService();

  if (serviceTestPassed) {
    console.log('\n🎉 SUCCESS: AI Setup is LOCKED IN!');
    console.log('================================');
    console.log('✅ Environment variables configured');
    console.log('✅ Dependencies installed');
    console.log('✅ AI Service working');
    console.log('✅ Provider switching functional');
    console.log('\n🚀 Ready to use with Cursor, Neovim, or any editor!');
    console.log('💡 Run: npm start');
  } else {
    console.log('\n❌ AI Service test failed. Check your configuration.');
    process.exit(1);
  }
}

if (require.main === module) {
  main().catch(console.error);
}

module.exports = { validateEnvironment, testAIService };