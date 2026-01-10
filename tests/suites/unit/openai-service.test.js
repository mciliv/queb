// Mock OpenAI at the module level
jest.mock('openai', () => {
  return {
    default: jest.fn().mockImplementation(() => ({
      chat: {
        completions: {
          create: jest.fn()
        }
      },
      completions: {
        create: jest.fn()
      }
    }))
  };
});

const AIService = require('../../../src/server/services/AIService');

describe('AIService LLM Call Coverage', () => {
  let aiService;
  let mockOpenAIClient;

  beforeEach(() => {
    // Get the mocked OpenAI client
    const OpenAI = require('openai').default;
    mockOpenAIClient = OpenAI.mock.results[0]?.value;

    if (!mockOpenAIClient) {
      // If no instance was created yet, create one
      OpenAI.mockClear();
      mockOpenAIClient = {
        chat: {
          completions: {
            create: jest.fn()
          }
        },
        completions: {
          create: jest.fn()
        }
      };
      OpenAI.mockImplementation(() => mockOpenAIClient);
    }

    // Clear all mocks
    jest.clearAllMocks();

    openAIService = new OpenAIService({
      apiKey: 'test-api-key',
      model: 'gpt-3.5-turbo',
      timeout: 5000,
      maxRetries: 2
    });
  });

  afterEach(() => {
    jest.clearAllMocks();
  });

  describe('Initialization', () => {
    test('should initialize with valid config', () => {
      expect(openAIService.client).toBeDefined();
      expect(openAIService.config.apiKey).toBe('test-api-key');
      expect(openAIService.config.model).toBe('gpt-3.5-turbo');
    });

    test('should throw error without API key', () => {
      expect(() => {
        new OpenAIService({ apiKey: null });
      }).toThrow('OpenAI API key is required');
    });
  });

  describe('callAPI - Chat Completions', () => {
    test('should successfully call chat completions API', async () => {
      const mockResponse = {
        choices: [{
          message: {
            content: '{"object": "water", "chemicals": [{"name": "Water", "smiles": "O"}]}'
          }
        }]
      };

      mockOpenAIClient.chat.completions.create.mockResolvedValue(mockResponse);

      const params = {
        messages: [{ role: 'user', content: 'Analyze water' }]
      };

      const result = await openAIService.callAPI(params);

      expect(mockOpenAIClient.chat.completions.create).toHaveBeenCalledWith(
        expect.objectContaining({
          model: 'gpt-3.5-turbo',
          messages: params.messages
        }),
        expect.any(Object)
      );

      expect(result).toEqual({
        object: 'water',
        chemicals: [{ name: 'Water', smiles: 'O' }]
      });
    });

    test('should use completion API when configured', async () => {
      const serviceWithCompletion = new OpenAIService({
        apiKey: 'test-api-key',
        model: 'gpt-3.5-turbo',
        useCompletionAPI: true
      });

      const mockResponse = {
        choices: [{
          text: '{"object": "water", "chemicals": [{"name": "Water", "smiles": "O"}]}'
        }]
      };

      mockOpenAIClient.completions.create.mockResolvedValue(mockResponse);

      const params = {
        messages: [{ role: 'user', content: 'Analyze water' }]
      };

      const result = await serviceWithCompletion.callAPI(params);

      expect(mockOpenAIClient.completions.create).toHaveBeenCalledWith(
        expect.objectContaining({
          model: 'text-davinci-003',
          prompt: expect.stringContaining('Analyze water')
        }),
        expect.any(Object)
      );

      expect(result).toEqual({
        object: 'water',
        chemicals: [{ name: 'Water', smiles: 'O' }]
      });
    });

    test('should handle API timeout', async () => {
      mockOpenAIClient.chat.completions.create.mockRejectedValue(
        new Error('Request timed out')
      );

      const params = {
        messages: [{ role: 'user', content: 'Analyze water' }]
      };

      await expect(openAIService.callAPI(params)).rejects.toThrow('Request timed out');
    });

    test('should handle rate limiting', async () => {
      const rateLimitError = new Error('Rate limit exceeded');
      rateLimitError.status = 429;
      rateLimitError.code = 'rate_limit_exceeded';

      mockOpenAIClient.chat.completions.create.mockRejectedValue(rateLimitError);

      const params = {
        messages: [{ role: 'user', content: 'Analyze water' }]
      };

      await expect(openAIService.callAPI(params)).rejects.toThrow('Rate limit exceeded');
    });

    test('should handle authentication errors', async () => {
      const authError = new Error('Invalid API key');
      authError.status = 401;
      authError.code = 'invalid_api_key';

      mockOpenAIClient.chat.completions.create.mockRejectedValue(authError);

      const params = {
        messages: [{ role: 'user', content: 'Analyze water' }]
      };

      await expect(openAIService.callAPI(params)).rejects.toThrow('Invalid API key');
    });

    test('should handle malformed JSON responses by returning raw content', async () => {
      const mockResponse = {
        choices: [{
          message: {
            content: '{invalid json response}'
          }
        }]
      };

      mockOpenAIClient.chat.completions.create.mockResolvedValue(mockResponse);

      const params = {
        messages: [{ role: 'user', content: 'Analyze water' }]
      };

      const result = await openAIService.callAPI(params);

      // The service returns raw content when JSON parsing fails
      expect(result).toBe('{invalid json response}');
    });

    test('should handle empty responses', async () => {
      const mockResponse = {
        choices: [{
          message: {
            content: ''
          }
        }]
      };

      mockOpenAIClient.chat.completions.create.mockResolvedValue(mockResponse);

      const params = {
        messages: [{ role: 'user', content: 'Analyze water' }]
      };

      await expect(openAIService.callAPI(params)).rejects.toThrow('Empty response content from AI');
    });

    test('should handle responses with no choices', async () => {
      const mockResponse = {
        choices: []
      };

      mockOpenAIClient.chat.completions.create.mockResolvedValue(mockResponse);

      const params = {
        messages: [{ role: 'user', content: 'Analyze water' }]
      };

      await expect(openAIService.callAPI(params)).rejects.toThrow('Invalid AI response structure');
    });

    test('should handle responses with null content', async () => {
      const mockResponse = {
        choices: [{
          message: {
            content: null
          }
        }]
      };

      mockOpenAIClient.chat.completions.create.mockResolvedValue(mockResponse);

      const params = {
        messages: [{ role: 'user', content: 'Analyze water' }]
      };

      await expect(openAIService.callAPI(params)).rejects.toThrow('Empty response content from AI');
    });

    test('should handle responses with undefined content', async () => {
      const mockResponse = {
        choices: [{
          message: {}
        }]
      };

      mockOpenAIClient.chat.completions.create.mockResolvedValue(mockResponse);

      const params = {
        messages: [{ role: 'user', content: 'Analyze water' }]
      };

      await expect(openAIService.callAPI(params)).rejects.toThrow('Empty response content from AI');
    });
  });

  describe('callAPI - Different Models', () => {
    test('should work with GPT-4 model', async () => {
      const gpt4Service = new OpenAIService({
        apiKey: 'test-api-key',
        model: 'gpt-4'
      });

      const mockResponse = {
        choices: [{
          message: {
            content: '{"result": "success"}'
          }
        }]
      };

      mockOpenAIClient.chat.completions.create.mockResolvedValue(mockResponse);

      const params = {
        messages: [{ role: 'user', content: 'Test' }]
      };

      await gpt4Service.callAPI(params);

      expect(mockOpenAIClient.chat.completions.create).toHaveBeenCalledWith(
        expect.objectContaining({ model: 'gpt-4' }),
        expect.any(Object)
      );
    });

    test('should work with GPT-5 model', async () => {
      const gpt5Service = new OpenAIService({
        apiKey: 'test-api-key',
        model: 'gpt-5'
      });

      const mockResponse = {
        choices: [{
          message: {
            content: '{"result": "success"}'
          }
        }]
      };

      mockOpenAIClient.chat.completions.create.mockResolvedValue(mockResponse);

      const params = {
        messages: [{ role: 'user', content: 'Test' }]
      };

      await gpt5Service.callAPI(params);

      expect(mockOpenAIClient.chat.completions.create).toHaveBeenCalledWith(
        expect.objectContaining({ model: 'gpt-5' }),
        expect.any(Object)
      );
    });
  });

  describe('Retry Logic', () => {
    test('should configure OpenAI client with max retries', () => {
      // The retry logic is handled by the OpenAI client itself
      // We can verify that the client is configured with the correct maxRetries
      const serviceWithRetries = new OpenAIService({
        apiKey: 'test-api-key',
        maxRetries: 5
      });

      // The OpenAI constructor should be called with maxRetries: 5
      const OpenAI = require('openai').default;
      expect(OpenAI).toHaveBeenCalledWith(
        expect.objectContaining({
          apiKey: 'test-api-key',
          timeout: expect.any(Number),
          maxRetries: 5
        })
      );
    });

    test('should use default max retries when not specified', () => {
      const serviceWithDefaults = new OpenAIService({
        apiKey: 'test-api-key'
      });

      // Should use the default maxRetries from config
      expect(serviceWithDefaults.config.maxRetries).toBe(2);
    });

    test('should handle errors from OpenAI client (retries handled internally)', async () => {
      const serverError = new Error('Internal server error');
      serverError.status = 500;

      mockOpenAIClient.chat.completions.create.mockRejectedValue(serverError);

      const params = {
        messages: [{ role: 'user', content: 'Test' }]
      };

      // The OpenAI client handles retries internally, so we just expect the final error
      await expect(openAIService.callAPI(params)).rejects.toThrow('Internal server error');
    });
  });

  describe('Health Check', () => {
    test('should return healthy status when API is accessible', async () => {
      mockOpenAIClient.chat.completions.create.mockResolvedValue({
        choices: [{
          message: {
            content: '{"status": "healthy"}'
          }
        }]
      });

      const health = await openAIService.healthCheck();

      expect(health).toEqual({
        status: 'healthy',
        response: {
          status: 'healthy'
        }
      });
    });

    test('should return unhealthy status when API is not accessible', async () => {
      mockOpenAIClient.chat.completions.create.mockRejectedValue(
        new Error('Connection failed')
      );

      const health = await openAIService.healthCheck();

      expect(health).toEqual({
        status: 'unhealthy',
        error: 'Connection failed'
      });
    });

    test('should return unhealthy status when API key is missing', async () => {
      expect(() => {
        new OpenAIService({
          apiKey: null
        });
      }).toThrow('OpenAI API key is required');
    });
  });

  describe('Response Parsing Edge Cases', () => {
    test('should handle deeply nested JSON', async () => {
      const complexResponse = {
        object: 'coffee',
        chemicals: [
          {
            name: 'Caffeine',
            smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
            properties: {
              molecular_weight: 194.19,
              formula: 'C8H10N4O2'
            }
          }
        ],
        metadata: {
          confidence: 0.95,
          source: 'PubChem'
        }
      };

      const mockResponse = {
        choices: [{
          message: {
            content: JSON.stringify(complexResponse)
          }
        }]
      };

      mockOpenAIClient.chat.completions.create.mockResolvedValue(mockResponse);

      const params = {
        messages: [{ role: 'user', content: 'Analyze coffee' }]
      };

      const result = await openAIService.callAPI(params);

      expect(result).toEqual(complexResponse);
    });

    test('should handle array responses', async () => {
      const arrayResponse = [
        { name: 'Water', smiles: 'O' },
        { name: 'Carbon Dioxide', smiles: 'O=C=O' }
      ];

      const mockResponse = {
        choices: [{
          message: {
            content: JSON.stringify(arrayResponse)
          }
        }]
      };

      mockOpenAIClient.chat.completions.create.mockResolvedValue(mockResponse);

      const params = {
        messages: [{ role: 'user', content: 'List chemicals' }]
      };

      const result = await openAIService.callAPI(params);

      expect(result).toEqual(arrayResponse);
    });

    test('should handle string responses', async () => {
      const stringResponse = 'This is a plain text response';

      const mockResponse = {
        choices: [{
          message: {
            content: stringResponse
          }
        }]
      };

      mockOpenAIClient.chat.completions.create.mockResolvedValue(mockResponse);

      const params = {
        messages: [{ role: 'user', content: 'Give me information' }]
      };

      const result = await openAIService.callAPI(params);

      expect(result).toBe(stringResponse);
    });

    test('should handle number responses', async () => {
      const numberResponse = 42;

      const mockResponse = {
        choices: [{
          message: {
            content: JSON.stringify(numberResponse)
          }
        }]
      };

      mockOpenAIClient.chat.completions.create.mockResolvedValue(mockResponse);

      const params = {
        messages: [{ role: 'user', content: 'What is the answer?' }]
      };

      const result = await openAIService.callAPI(params);

      expect(result).toBe(numberResponse);
    });

    test('should handle boolean responses', async () => {
      const booleanResponse = true;

      const mockResponse = {
        choices: [{
          message: {
            content: JSON.stringify(booleanResponse)
          }
        }]
      };

      mockOpenAIClient.chat.completions.create.mockResolvedValue(mockResponse);

      const params = {
        messages: [{ role: 'user', content: 'Is this correct?' }]
      };

      const result = await openAIService.callAPI(params);

      expect(result).toBe(booleanResponse);
    });
  });

  describe('Development Logging', () => {
    beforeEach(() => {
      process.env.NODE_ENV = 'development';
    });

    afterEach(() => {
      delete process.env.NODE_ENV;
    });

    test('should log request details in development', async () => {
      const consoleSpy = jest.spyOn(console, 'log').mockImplementation();

      const mockResponse = {
        choices: [{
          message: {
            content: '{"result": "success"}'
          }
        }]
      };

      mockOpenAIClient.chat.completions.create.mockResolvedValue(mockResponse);

      const params = {
        messages: [{ role: 'user', content: 'Test' }],
        max_tokens: 100
      };

      await openAIService.callAPI(params);

      expect(consoleSpy).toHaveBeenCalledWith(
        '[OpenAIService] Chat Completion request:',
        expect.objectContaining({
          model: 'gpt-3.5-turbo',
          messagesCount: 1,
          maxTokens: 100
        })
      );

      expect(consoleSpy).toHaveBeenCalledWith('[OpenAIService] Chat Completion response received');

      consoleSpy.mockRestore();
    });

    test('should not log in production', async () => {
      process.env.NODE_ENV = 'production';
      const consoleSpy = jest.spyOn(console, 'log').mockImplementation();

      const mockResponse = {
        choices: [{
          message: {
            content: '{"result": "success"}'
          }
        }]
      };

      mockOpenAIClient.chat.completions.create.mockResolvedValue(mockResponse);

      const params = {
        messages: [{ role: 'user', content: 'Test' }]
      };

      await openAIService.callAPI(params);

      expect(consoleSpy).not.toHaveBeenCalled();

      consoleSpy.mockRestore();
    });
  });

  describe('Environment Variables', () => {
    test('should use completion API when USE_COMPLETION_API is set', async () => {
      process.env.USE_COMPLETION_API = 'true';

      const serviceWithEnv = new OpenAIService({
        apiKey: 'test-api-key',
        model: 'gpt-3.5-turbo'
      });

      const mockResponse = {
        choices: [{
          text: '{"result": "success"}'
        }]
      };

      mockOpenAIClient.completions.create.mockResolvedValue(mockResponse);

      const params = {
        messages: [{ role: 'user', content: 'Test' }]
      };

      const result = await serviceWithEnv.callAPI(params);

      expect(mockOpenAIClient.completions.create).toHaveBeenCalled();
      expect(result).toEqual({ result: 'success' });

      delete process.env.USE_COMPLETION_API;
    });

    test('should prefer options over environment variable', async () => {
      process.env.USE_COMPLETION_API = 'true';

      const mockResponse = {
        choices: [{
          message: {
            content: '{"result": "success"}'
          }
        }]
      };

      mockOpenAIClient.chat.completions.create.mockResolvedValue(mockResponse);

      const params = {
        messages: [{ role: 'user', content: 'Test' }]
      };

      // Force chat completions despite env var
      const result = await openAIService.callAPI(params, { useCompletionAPI: false });

      expect(mockOpenAIClient.chat.completions.create).toHaveBeenCalled();
      expect(result).toEqual({ result: 'success' });

      delete process.env.USE_COMPLETION_API;
    });
  });
});
