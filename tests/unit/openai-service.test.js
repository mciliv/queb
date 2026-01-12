// Kept under the legacy filename for continuity, but updated to the current architecture.
// This is a true unit test for `AIService` using mocked `@ai-sdk` + `ai.generateText`.

jest.mock('ai', () => ({
  generateText: jest.fn(),
}));

jest.mock('@ai-sdk/openai', () => ({
  openai: jest.fn(() => ({ provider: 'openai', modelId: 'mock-openai-model' })),
}));

jest.mock('@ai-sdk/xai', () => ({
  xai: jest.fn(() => ({ provider: 'xai', modelId: 'mock-xai-model' })),
}));

const { generateText } = require('ai');
const AIService = require('../../src/server/services/AIService');

describe('AIService (legacy openai-service.test.js) - unit', () => {
  const originalEnv = process.env;

  beforeEach(() => {
    jest.clearAllMocks();
    process.env = { ...originalEnv };
    process.env.NODE_ENV = 'test';
    process.env.AI_PROVIDER = 'openai';
    process.env.OPENAI_API_KEY = 'test-openai-key';
    process.env.OPENAI_MODEL = 'latest';
  });

  afterAll(() => {
    process.env = originalEnv;
  });

  test('throws if OPENAI_API_KEY is missing', () => {
    delete process.env.OPENAI_API_KEY;
    expect(() => new AIService()).toThrow('OPENAI API key is required');
  });

  test('parses JSON response text', async () => {
    generateText.mockResolvedValue({
      text: '{"object":"water","chemicals":[{"name":"Water","smiles":"O"}]}',
      finishReason: 'stop',
      usage: { promptTokens: 1, completionTokens: 1, totalTokens: 2 },
    });

    const svc = new AIService();
    const result = await svc.callAPI({ messages: [{ role: 'user', content: 'Analyze water' }] });

    expect(result).toEqual({ object: 'water', chemicals: [{ name: 'Water', smiles: 'O' }] });
  });

  test('returns structured object when response is not JSON', async () => {
    generateText.mockResolvedValue({
      text: 'plain text response',
      finishReason: 'stop',
      usage: { promptTokens: 1, completionTokens: 1, totalTokens: 2 },
    });

    const svc = new AIService();
    const result = await svc.callAPI({ messages: [{ role: 'user', content: 'hi' }] });

    expect(result).toEqual(
      expect.objectContaining({
        content: 'plain text response',
        role: 'assistant',
      }),
    );
  });
});

