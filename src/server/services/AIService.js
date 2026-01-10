/**
 * Unified AI Service for multiple providers (OpenAI, xAI, etc.)
 * Uses official SDKs instead of hard-coded URLs
 * Always loads configuration from environment variables
 */
class AIService {
  constructor() {
    // Load configuration from environment variables (.env file)
    const isTest = process.env.NODE_ENV === 'test';

    this.config = {
      provider: process.env.AI_PROVIDER || 'openai'
    };

    // Provider-specific configurations loaded from .env
    // Test environment always uses mock-model regardless of .env settings
    this.providerConfigs = {
      openai: {
        apiKey: process.env.OPENAI_API_KEY,
        model: isTest ? 'mock-model' : (process.env.OPENAI_MODEL || 'latest'),
        timeout: isTest ? 1000 : 30000
      },
      xai: {
        apiKey: process.env.XAI_API_KEY,
        model: isTest ? 'mock-model' : (process.env.XAI_MODEL || 'latest'),
        timeout: isTest ? 1000 : 30000
      }
    };

    this.client = null;
    this.model = null;
    this.initializeClient();
  }

  /**
   * Get default model for a provider
   */
  _getDefaultModel() {
    const defaults = {
      openai: 'gpt-3.5-turbo',
      xai: 'grok-3'
    };
    return defaults[this.config.provider] || 'gpt-3.5-turbo';
  }

  /**
   * Get default base URL for a provider
   */
  _getDefaultBaseURL() {
    const defaults = {
      openai: 'https://api.openai.com/v1',
      xai: 'https://api.x.ai/v1'
    };
    return defaults[this.config.provider] || 'https://api.openai.com/v1';
  }

  /**
   * Get current provider config
   */
  _getCurrentProviderConfig() {
    return this.providerConfigs[this.config.provider] || this.providerConfigs.openai;
  }

  /**
   * Initialize the client and model for the current provider using official SDKs
   */
  initializeClient() {
    const providerConfig = this._getCurrentProviderConfig();

    if (!providerConfig.apiKey) {
      throw new Error(`${this.config.provider.toUpperCase()} API key is required`);
    }

    // Use official SDKs instead of hard-coded URLs
    switch (this.config.provider) {
      case 'openai':
        this.initializeOpenAIClient(providerConfig);
        break;
      case 'xai':
        this.initializeXAIClient(providerConfig);
        break;
      default:
        throw new Error(`Unsupported provider: ${this.config.provider}`);
    }
  }

  /**
   * Initialize OpenAI client using official SDK
   */
  initializeOpenAIClient(config) {
    const { openai } = require('@ai-sdk/openai');

    this.model = openai(config.model, {
      apiKey: config.apiKey
    });

    // Keep a reference to the raw client for advanced operations
    const OpenAI = require('openai').default || require('openai');
    this.client = new OpenAI({
      apiKey: config.apiKey,
      timeout: config.timeout,
      maxRetries: 2
    });
  }

  /**
   * Initialize xAI client using official SDK
   */
  initializeXAIClient(config) {
    const { xai } = require('@ai-sdk/xai');

    this.model = xai(config.model, {
      apiKey: config.apiKey
    });

    // Keep a reference to the raw client for advanced operations
    // xAI SDK provides the model interface, but we may need raw client for some operations
    const OpenAI = require('openai').default || require('openai');
    this.client = new OpenAI({
      apiKey: config.apiKey,
      baseURL: 'https://api.x.ai/v1', // Official xAI endpoint
      timeout: config.timeout,
      maxRetries: 2
    });
  }

  /**
   * Generic method to call AI API using Responses API
   */
  async callAPI(params, options = {}) {
    // Use the Responses API
    const result = await this._callResponsesAPI(params);
    return this._parseResponse(result);
  }

  /**
   * Call Responses API
   */
  async _callResponsesAPI(params) {
    // Handle mock model for testing
    if (this._getCurrentProviderConfig().model === 'mock-model') {
      return this._mockResponse(params);
    }

    const providerConfig = this._getCurrentProviderConfig();

    // Convert chat messages to Responses API format
    let requestParams;

    if (params.messages && Array.isArray(params.messages)) {
      // Convert messages array to input string
      const input = params.messages.map(msg => {
        const role = msg.role || 'user';
        const content = msg.content || '';
        return `${role.toUpperCase()}: ${content}`;
      }).join('\n\n');

      requestParams = {
        model: providerConfig.model,
        input: input,
        ...params
      };

      // Remove messages since we converted to input
      delete requestParams.messages;
    } else {
      // Use input directly if provided
      requestParams = {
        model: providerConfig.model,
        input: params.input || params.prompt || '',
        ...params
      };
    }

    if (process.env.NODE_ENV === 'development') {
      console.log(`[AIService:${this.config.provider}] Responses API request:`, {
        model: requestParams.model,
        inputLength: requestParams.input?.length,
        maxTokens: requestParams.max_tokens
      });
    }

    try {
      const result = await this.client.responses.create(
        requestParams,
        { timeout: providerConfig.timeout || 30000 }
      );

      if (process.env.NODE_ENV === 'development') {
        console.log(`[AIService:${this.config.provider}] Responses API response received`);
      }

      return result;
    } catch (error) {
      console.error(`[AIService:${this.config.provider}] Responses API error:`, error.message);
      throw error;
    }
  }


  /**
   * Parse response from Responses API
   */
  _parseResponse(response) {
    try {
      // Check for Responses API format first (output_text)
      if (response?.output_text) {
        const content = response.output_text;

        if (!content.trim()) {
          throw new Error('Empty response content from AI');
        }

        // Try to parse as JSON first
        try {
          return JSON.parse(content);
        } catch {
          // Return structured response if not JSON
          return {
            content: content,
            role: 'assistant',
            finish_reason: response.finish_reason || 'stop',
            usage: response.usage,
            id: response.id,
            model: response.model
          };
        }
      }

      // Fallback to chat completions format for backward compatibility
      if (response?.choices?.[0]) {
        const choice = response.choices[0];
        const content = choice.message?.content || choice.text || '';

        if (!content.trim()) {
          throw new Error('Empty response content from AI');
        }

        // Try to parse as JSON first
        try {
          return JSON.parse(content);
        } catch {
          // Return structured response if not JSON
          return {
            content: content,
            role: choice.message?.role || 'assistant',
            finish_reason: choice.finish_reason,
            usage: response.usage,
            id: response.id,
            model: response.model
          };
        }
      }

      throw new Error('Invalid AI response structure - no output_text or choices found');
    } catch (error) {
      if (process.env.NODE_ENV === 'development') {
        console.error(`[AIService:${this.config.provider}] Failed to parse response:`, {
          error: error.message,
          responseStructure: {
            hasOutputText: !!response?.output_text,
            hasChoices: !!response?.choices,
            choiceCount: response?.choices?.length,
            firstChoice: response?.choices?.[0] ? {
              hasMessage: !!response.choices[0].message,
              hasText: !!response.choices[0].text,
              hasContent: !!(response.choices[0].message?.content || response.choices[0].text)
            } : null
          }
        });
      }
      throw error;
    }
  }


  /**
   * Mock response for testing (Responses API format)
   */
  _mockResponse(params) {
    return {
      id: 'mock-response-id',
      object: 'response',
      created: Date.now(),
      model: 'mock-model',
      output_text: 'This is a mock response for testing purposes.',
      finish_reason: 'stop',
      usage: {
        prompt_tokens: 10,
        completion_tokens: 20,
        total_tokens: 30
      }
    };
  }

  /**
   * Get current provider name
   */
  getProvider() {
    return this.config.provider;
  }

  /**
   * Get current model
   */
  getModel() {
    return this._getCurrentProviderConfig().model;
  }

  /**
   * Switch provider at runtime (reloads config from .env)
   */
  switchProvider(newProvider) {
    if (!this.providerConfigs[newProvider]) {
      throw new Error(`Unsupported provider: ${newProvider}. Supported: ${Object.keys(this.providerConfigs).join(', ')}`);
    }

    // Switch provider and reload all config from .env
    const isTest = process.env.NODE_ENV === 'test';

    this.config.provider = newProvider;
    this.providerConfigs = {
      openai: {
        apiKey: process.env.OPENAI_API_KEY,
        model: process.env.OPENAI_MODEL || (isTest ? 'mock-model' : 'latest'),
        timeout: isTest ? 1000 : 30000
      },
      xai: {
        apiKey: process.env.XAI_API_KEY,
        model: process.env.XAI_MODEL || 'latest',
        timeout: isTest ? 1000 : 30000
      }
    };

    // Re-initialize client with new provider
    this.initializeClient();
  }

  /**
   * Health check to validate API connectivity and contract
   */
  async healthCheck() {
    try {
      const response = await this.callAPI({
        messages: [{ role: 'user', content: 'Hello' }],
        max_tokens: 5
      });
      return { status: 'healthy', response };
    } catch (error) {
      return { status: 'unhealthy', error: error.message };
    }
  }

  /**
   * Get available providers
   */
  static getAvailableProviders() {
    return ['openai', 'xai'];
  }

  /**
   * Get provider-specific defaults
   */
  static getProviderDefaults(provider) {
    const defaults = {
      openai: {
        baseURL: 'https://api.openai.com/v1',
        model: 'gpt-3.5-turbo'
      },
      xai: {
        baseURL: 'https://api.x.ai/v1',
        model: 'grok-3'
      }
    };
    return defaults[provider] || defaults.openai;
  }
}

module.exports = AIService;