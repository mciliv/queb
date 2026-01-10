#!/usr/bin/env node

/**
 * Test script to compare OpenAI Chat Completions vs Completion API
 * Run with: USE_COMPLETION_API=true node test-completion-api.js
 */

const OpenAI = require('openai').default || require('openai');

async function testAPIs() {
  const useCompletionAPI = process.env.USE_COMPLETION_API === 'true';
  const apiKey = process.env.OPENAI_API_KEY;

  if (!apiKey) {
    console.error('❌ OPENAI_API_KEY environment variable required');
    process.exit(1);
  }

  const client = new OpenAI({ apiKey });

  const testPrompt = `You are a chemical analysis expert. Analyze the given object and identify the chemical compounds present.

For each chemical compound, provide:
- name: Common chemical name
- smiles: SMILES notation if known

Return your analysis as a JSON object with this structure:
{
  "object": "description of the analyzed object",
  "chemicals": [
    {
      "name": "chemical name",
      "smiles": "SMILES notation"
    }
  ]
}

Object: coffee`;

  console.log(`🧪 Testing ${useCompletionAPI ? 'COMPLETION' : 'CHAT COMPLETION'} API`);
  console.log('═'.repeat(60));

  try {
    let response;
    const startTime = Date.now();

    if (useCompletionAPI) {
      // Use Completion API
      response = await client.completions.create({
        model: 'text-davinci-003', // or gpt-3.5-turbo-instruct
        prompt: testPrompt,
        max_tokens: 200,
        temperature: 1
      });
    } else {
      // Use Chat Completion API (default)
      response = await client.chat.completions.create({
        model: 'gpt-3.5-turbo',
        messages: [{ role: 'user', content: testPrompt }],
        max_tokens: 200,
        temperature: 1
      });
    }

    const duration = Date.now() - startTime;

    console.log(`✅ Response received in ${duration}ms`);
    console.log('Response structure:', {
      choices: response.choices?.length,
      usage: response.usage
    });

    // Extract content
    const content = useCompletionAPI
      ? response.choices[0]?.text
      : response.choices[0]?.message?.content;

    console.log('\n📄 Raw AI Response:');
    console.log(content);

    // Try to parse as JSON
    try {
      const parsed = JSON.parse(content);
      console.log('\n📋 Parsed JSON:');
      console.log(JSON.stringify(parsed, null, 2));
    } catch (parseError) {
      console.log('\n⚠️  Response is not valid JSON:', parseError.message);
    }

  } catch (error) {
    console.error('❌ API call failed:', error.message);
    if (error.response) {
      console.log('Error details:', error.response.data);
    }
  }
}

testAPIs();


