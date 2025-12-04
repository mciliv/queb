/**
 * Unit tests for AI service (OpenAI client)
 */

const { config, createOpenAIClient } = require('../../../src/core/services');
const { OpenAI } = require('openai');

describe('AI Service', () => {

  describe('Environment Check', () => {
    test('should show current environment state', () => {
      const debugInfo = config.getDebugInfo();
      console.log('Configuration debug:', JSON.stringify(debugInfo, null, 2));
      console.log('AI Config apiKey:', config.openai.apiKey ? `SET (length: ${config.openai.apiKey.length})` : 'NOT SET');
      console.log('process.env.OPENAI_API_KEY:', process.env.OPENAI_API_KEY ? `SET (length: ${process.env.OPENAI_API_KEY.length})` : 'NOT SET');
    });
  });

  describe('Configuration', () => {
    test('should load AI config with required fields', () => {
      expect(config.openai).toHaveProperty('apiKey');
      expect(config.openai).toHaveProperty('model');
      expect(config.openai).toHaveProperty('timeout');
    });

    test('should have valid model specified', () => {
      expect(config.openai.model).toBeTruthy();
      expect(typeof config.openai.model).toBe('string');
    });

    test('should have valid timeout', () => {
      expect(config.openai.timeout).toBeTruthy();
      expect(typeof config.openai.timeout).toBe('number');
    });
  });

  describe('Client Creation', () => {
    test('should create client successfully', () => {
      const client = createOpenAIClient();

      if (!OpenAI) {
        expect(client).toBeNull();
      } else {
        expect(client).toBeTruthy();
        expect(client).toBeInstanceOf(OpenAI);
      }
    });

    test('should create client with proper configuration', () => {
      const client = createOpenAIClient();

      if (client) {
        // Verify client has expected OpenAI methods
        expect(typeof client.chat?.completions?.create).toBe('function');
      }
    });
  });

  describe('Client Functionality', () => {
    let client;

    beforeEach(() => {
      client = createOpenAIClient();
    });

    test('should have chat completions API', () => {
      if (client) {
        expect(client.chat).toBeDefined();
        expect(client.chat.completions).toBeDefined();
        expect(typeof client.chat.completions.create).toBe('function');
      }
    });

    test('should make successful API call with valid key', async () => {
      if (!client || !process.env.OPENAI_API_KEY) {
        console.log('⏭️  Skipping live API test (no client or API key)');
        return;
      }

      try {
        const response = await client.chat.completions.create({
          model: config.openai.model,
          messages: [
            { role: 'system', content: 'You are a chemistry expert. Respond in one sentence.' },
            { role: 'user', content: 'What is the molecular formula for water?' }
          ],
          temperature: 0.2,
          max_tokens: 50
        });

        expect(response).toBeDefined();
        expect(response.choices).toBeDefined();
        expect(response.choices.length).toBeGreaterThan(0);
        expect(response.choices[0].message).toBeDefined();
        expect(response.choices[0].message.content).toBeTruthy();

        console.log('✅ AI response:', response.choices[0].message.content);
      } catch (error) {
        if (error.status === 401) {
          console.log('⏭️  Skipping live API test (invalid API key)');
        } else {
          throw error;
        }
      }
    }, 30000); // 30 second timeout for API call
  });
});