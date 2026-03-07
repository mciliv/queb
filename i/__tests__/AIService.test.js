const AIService = require('../AIService');

describe('AIService', () => {
  const originalEnv = process.env;

  beforeEach(() => {
    process.env = {
      ...originalEnv,
      NODE_ENV: 'test',
      AI_PROVIDER: 'openai',
      OPENAI_API_KEY: 'test-key',
      OPENAI_MODEL: ''
    };
  });

  afterEach(() => {
    process.env = originalEnv;
  });

  it('returns normalized content for mock responses', async () => {
    const service = new AIService();
    const result = await service.callAPI({
      input: 'Hello',
      max_tokens: 10
    });

    expect(result).toEqual(
      expect.objectContaining({
        content: 'This is a mock response for testing purposes.',
        role: 'assistant',
        finish_reason: 'stop',
        model: 'mock-model'
      })
    );
  });
});
