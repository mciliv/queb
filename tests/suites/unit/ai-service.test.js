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
const AIService = require('../../../src/server/services/AIService');

describe('AIService (@ai-sdk wrapper) - unit', () => {
  const originalEnv = process.env;

  beforeEach(() => {
    jest.clearAllMocks();
    process.env = { ...originalEnv };
    process.env.NODE_ENV = 'test';
    process.env.AI_PROVIDER = 'openai';
    process.env.OPENAI_API_KEY = 'test-openai-key';
  });

  afterAll(() => {
    process.env = originalEnv;
  });

  test('throws if provider API key is missing', () => {
    delete process.env.OPENAI_API_KEY;
    expect(() => new AIService()).toThrow('OPENAI API key is required');
  });

  test('parses JSON from @ai-sdk response text', async () => {
    process.env.OPENAI_MODEL = 'latest';

    generateText.mockResolvedValue({
      text: '{"object":"water","chemicals":[{"name":"Water","smiles":"O"}]}',
      usage: { promptTokens: 1, completionTokens: 1, totalTokens: 2 },
      finishReason: 'stop',
    });

    const svc = new AIService();
    const result = await svc.callAPI({ messages: [{ role: 'user', content: 'Analyze water' }] });

    expect(generateText).toHaveBeenCalledWith(
      expect.objectContaining({
        model: expect.any(Object),
        messages: [{ role: 'user', content: 'Analyze water' }],
      }),
    );

    expect(result).toEqual({
      object: 'water',
      chemicals: [{ name: 'Water', smiles: 'O' }],
    });
  });

  test('returns structured object when response is not JSON', async () => {
    process.env.OPENAI_MODEL = 'latest';

    generateText.mockResolvedValue({
      text: 'hello world',
      usage: { promptTokens: 1, completionTokens: 1, totalTokens: 2 },
      finishReason: 'stop',
    });

    const svc = new AIService();
    const result = await svc.callAPI({ messages: [{ role: 'user', content: 'hi' }] });

    expect(result).toEqual(
      expect.objectContaining({
        content: 'hello world',
        role: 'assistant',
      }),
    );
  });

  test('propagates SDK errors (e.g. timeout)', async () => {
    process.env.OPENAI_MODEL = 'latest';

    generateText.mockRejectedValue(new Error('Request timed out.'));
    const svc = new AIService();
    await expect(svc.callAPI({ messages: [{ role: 'user', content: 'x' }] })).rejects.toThrow(
      'Request timed out.',
    );
  });
});

