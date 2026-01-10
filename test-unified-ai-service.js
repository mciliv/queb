#!/usr/bin/env node

/**
 * Test script to demonstrate the unified AI service working with both OpenAI and xAI
 */

const UnifiedAIService = require('./src/server/services/UnifiedAIService');

async function testUnifiedAIService() {
  console.log('🧪 Testing Unified AI Service\n');

  // Test data
  const testParams = {
    messages: [
      { role: 'user', content: 'Hello, can you tell me a short joke?' }
    ],
    max_tokens: 50,
    temperature: 0.7
  };

  // Test OpenAI (mock mode for demo)
  console.log('1️⃣ Testing OpenAI Provider (mock mode):');
  try {
    const openAIService = new UnifiedAIService({
      provider: 'openai',
      model: 'mock-model', // Use mock for testing
      apiKey: 'test-key'
    });

    console.log(`   Provider: ${openAIService.getProvider()}`);
    console.log(`   Model: ${openAIService.getModel()}`);

    const openAIResult = await openAIService.callAPI(testParams);
    console.log(`   ✅ Response: ${openAIResult.content}`);
    console.log(`   Role: ${openAIResult.role}`);
    console.log(`   Model: ${openAIResult.model}\n`);
  } catch (error) {
    console.log(`   ❌ Error: ${error.message}\n`);
  }

  // Test xAI (mock mode for demo)
  console.log('2️⃣ Testing xAI Provider (mock mode):');
  try {
    const xaiService = new UnifiedAIService({
      provider: 'xai',
      model: 'mock-model', // Use mock for testing
      apiKey: 'test-key'
    });

    console.log(`   Provider: ${xaiService.getProvider()}`);
    console.log(`   Model: ${xaiService.getModel()}`);

    const xaiResult = await xaiService.callAPI(testParams);
    console.log(`   ✅ Response: ${xaiResult.content}`);
    console.log(`   Role: ${xaiResult.role}`);
    console.log(`   Model: ${xaiResult.model}\n`);
  } catch (error) {
    console.log(`   ❌ Error: ${error.message}\n`);
  }

  // Test provider switching
  console.log('3️⃣ Testing Provider Switching:');
  try {
    const service = new UnifiedAIService({
      provider: 'openai',
      model: 'mock-model',
      apiKey: 'test-key'
    });

    console.log(`   Initial provider: ${service.getProvider()}`);

    service.switchProvider('xai', {
      model: 'mock-model',
      apiKey: 'test-key'
    });

    console.log(`   Switched to provider: ${service.getProvider()}`);

    const switchedResult = await service.callAPI(testParams);
    console.log(`   ✅ Response after switch: ${switchedResult.content}\n`);
  } catch (error) {
    console.log(`   ❌ Error: ${error.message}\n`);
  }

  // Show available providers
  console.log('4️⃣ Available Providers:');
  console.log(`   ${UnifiedAIService.getAvailableProviders().join(', ')}\n`);

  // Show provider defaults
  console.log('5️⃣ Provider Defaults:');
  UnifiedAIService.getAvailableProviders().forEach(provider => {
    const defaults = UnifiedAIService.getProviderDefaults(provider);
    console.log(`   ${provider}: ${defaults.model} (${defaults.baseURL})`);
  });

  console.log('\n🎉 Unified AI Service test completed successfully!');
  console.log('\n💡 To use with real APIs, set these environment variables:');
  console.log('   OpenAI: OPENAI_API_KEY, XAI_API_KEY');
  console.log('   Provider: AI_PROVIDER=openai or AI_PROVIDER=xai');
}

// Run the test
if (require.main === module) {
  testUnifiedAIService().catch(console.error);
}

module.exports = { testUnifiedAIService };
