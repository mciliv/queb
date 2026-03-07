/**
 * AIService architectural constraints (kept close to code):
 * - **Config**: `.env` / `process.env` is the single source of truth → see constructor + `_createProviderConfigs()`
 * - **Execution**: use `@ai-sdk` (`generateText`) only → see `callAPI()` + `_callWithSDK()`
 * - **Conversion/Parsing**: internal ↔ SDK formats are centralized → see `_convertToSDKFormat()` + `_parseSDKResponse()`
 *
 * If you need to change behavior, change it at those specific choke points (don’t scatter config or raw SDK calls).
 */
class AIService {
  constructor() {
    /**
     * SINGLE SOURCE OF TRUTH (CONFIG)
     * - All AI config must come from `.env` → `process.env`
     * - No constructor options, no hardcoded keys/models for runtime behavior
     */

    this.config = {
      provider: process.env.AI_PROVIDER || 'openai'
    };

    // Provider-specific configurations loaded from .env.
    // Test env uses `mock-model` to prevent network calls.
    this.providerConfigs = this._createProviderConfigs();

    this.model = null;

    /**
     * AGENTIC RESILIENCE LAYER
     * - ModelManager handles 429 failover
     * - AgentRegistry handles task-specialized model priority
     */
    const ModelManager = require('./lib/ModelManager');
    const AgentRegistry = require('./lib/AgentRegistry');
    this.modelManager = new ModelManager(this);
    this.agentRegistry = new AgentRegistry(this.modelManager);

    /**
     * SDK ABSTRACTION PURITY
     * - Initialize only an `@ai-sdk` model adapter
     * - Do NOT initialize or hold raw provider clients here
     */
    this.initializeModel();
  }

  /**
   * Run specialized agent task
   */
  async runAgentTask(agent, prompt, params = {}) {
    return await this.agentRegistry.runTask(agent, prompt, params);
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
   * CONFIGURATION FACTORY - SINGLE SOURCE OF TRUTH
   *
   * This method is the EXCLUSIVE source for all AI provider configuration.
   * Every config decision flows through here.
   *
   * REQUIRED BEHAVIORS:
   * - Read ONLY from process.env
   - Return consistent config objects
   - Handle test/prod differences here only
   *
   * WHY CENTRALIZED:
   * - Prevents config drift across the codebase
   - Enables atomic config changes
   - Maintains provider abstraction
   - Leverages SDK reliability features (auto model versioning, rate limiting, retries)
   *
   * USAGE: Call this method whenever provider config is needed
   */
  _createProviderConfigs() {
    return {
      openai: this._createProviderConfig('OPENAI'),
      xai: this._createProviderConfig('XAI')
    };
  }

  /**
   * PROVIDER CONFIGURATION HELPER - DRY CONFIGURATION
   *
   * Creates standardized configuration for any AI provider.
   * Eliminates duplication and ensures consistency.
   *
   * @param {string} prefix - Environment variable prefix ('OPENAI' or 'XAI')
   * @returns {object} Provider configuration object
   */
  _createProviderConfig(prefix) {
    /**
     * CONFIG FACTORY
     * - Centralize env var reads here to avoid config drift
     * - Keep test/prod differences here (not scattered across call sites)
     */
    const isTest = process.env.NODE_ENV === 'test';
    const envKey = `${prefix}_API_KEY`;
    const envModel = `${prefix}_MODEL`;

    return {
      apiKey: process.env[envKey],
      // SDK automatically resolves 'latest' to newest available model
      model: process.env[envModel] || (isTest ? 'mock-model' : 'latest'),
      timeout: isTest ? 1000 : 30000
    };
  }

  /**
   * PROVIDER MODEL FACTORY - MODULAR INITIALIZATION
   *
   * Registry of provider-specific model initialization functions.
   * Each provider has its own initialization method for clean separation.
   */
  _getProviderModelFactory() {
    return {
      openai: (config) => this._initializeOpenAIModel(config),
      xai: (config) => this._initializeXAIModel(config),
      google: (config) => this._initializeGoogleModel(config)
    };
  }

  /**
   * OPENAI MODEL INITIALIZATION - PROVIDER SPECIFIC
   *
   * Initializes OpenAI model with @ai-sdk/openai adapter.
   * Isolated for clean modular design and easy testing.
   *
   * @param {object} config - Provider configuration with apiKey and model
   * @returns {@ai-sdk model} Initialized OpenAI model
   */
  _initializeOpenAIModel(config) {
    const { openai } = require('@ai-sdk/openai');
    return openai(config.model, {
      apiKey: config.apiKey
    });
  }

  /**
   * XAI MODEL INITIALIZATION - PROVIDER SPECIFIC
   *
   * Initializes xAI model with @ai-sdk/xai adapter.
   * Isolated for clean modular design and easy testing.
   *
   * @param {object} config - Provider configuration with apiKey and model
   * @returns {@ai-sdk model} Initialized xAI model
   */
  _initializeXAIModel(config) {
    const { xai } = require('@ai-sdk/xai');
    return xai(config.model, {
      apiKey: config.apiKey
    });
  }

  /**
   * GOOGLE MODEL INITIALIZATION - PROVIDER SPECIFIC
   *
   * Initializes Google model with @ai-sdk/google adapter.
   * Isolated for clean modular design and easy testing.
   *
   * @param {object} config - Provider configuration with apiKey and model
   * @returns {@ai-sdk model} Initialized Google model
   */
  _initializeGoogleModel(config) {
    const { google } = require('@ai-sdk/google');
    return google(config.model, {
      apiKey: config.apiKey
    });
  }

  /**
   * MODEL INITIALIZATION - SDK ABSTRACTION LAYER
   *
   * Creates @ai-sdk model instances for unified AI interactions.
   * This is the ONLY initialization needed - no raw clients required.
   *
   * BENEFITS:
   * - Automatic provider abstraction
   * - Consistent API across all providers
   * - Future-proof against SDK changes
   * - Modular provider handling
   *
   * RESULT: this.model ready for generateText() calls
   */
  initializeModel() {
    /**
     * MODEL INIT (ABSTRACTION)
     * - Provider switch is allowed, but the runtime object must still be an `@ai-sdk` model
     * - The only provider-specific code allowed here is choosing the adapter factory.
     */
    const providerConfig = this._getCurrentProviderConfig();

    if (!providerConfig.apiKey) {
      throw new Error(`${this.config.provider.toUpperCase()} API key is required`);
    }

    // MODULAR PROVIDER INITIALIZATION
    // Use factory pattern for clean separation of provider-specific logic
    const providerFactories = this._getProviderModelFactory();

    if (!providerFactories[this.config.provider]) {
      throw new Error(`Unsupported provider: ${this.config.provider}`);
    }

    this.model = providerFactories[this.config.provider](providerConfig);
  }

  /**
   * UNIFIED AI INTERFACE - EXTERNAL API
   *
   * This is the ONLY public method for AI interactions.
   * All complexity is abstracted away - callers get consistent results.
   *
   * INPUT: Standard params (messages, input, temperature, etc.)
   * OUTPUT: Normalized response (content, usage, metadata)
   *
   * WHY UNIFIED:
   * - Same interface regardless of provider
   * - Automatic format conversion
   * - Consistent error handling
   * - Provider switching transparency
   */
  async callAPI(params, options = {}) {
    try {
      const result = await this._callWithSDK(params);
      return this._parseSDKResponse(result);
    } catch (error) {
      throw error;
    }
  }

  /**
   * SDK EXECUTION LAYER - CORE AI INTERACTION
   *
   * This method contains the actual @ai-sdk generateText() call.
   * It's the bridge between internal format and SDK expectations.
   *
   * FLOW:
   * 1. Check for test mocks
   * 2. Convert internal params to SDK format
   * 3. Call generateText()
   * 4. Return unified response
   *
   * WHY ISOLATED:
   * - Contains all SDK-specific logic
   * - Easy to test and mock
   * - Single point for SDK updates
   * - Leverages SDK reliability: automatic retries, rate limiting, error recovery
   */
  async _callWithSDK(params) {
    const { generateText } = require('ai');

    // Handle mock model for testing
    if (this._getCurrentProviderConfig().model === 'mock-model') {
      return this._getMockResponse(params);
    }

    /**
     * SDK RELIABILITY MAXIMIZATION
     * - Prefer `generateText({ model, ... })` so retries/compat are handled by the SDK stack
     * - Keep parameter mapping centralized in `_convertToSDKFormat()`
     */
    // Convert params to @ai-sdk format
    const generateParams = this._convertToSDKFormat(params);

    if (process.env.NODE_ENV === 'development') {
      console.log(`[AIService:${this.config.provider}] @ai-sdk generateText request:`, {
        model: this.model.modelId || this.model,
        inputLength: generateParams.prompt?.length || generateParams.messages?.length,
        maxTokens: generateParams.maxTokens
      });
    }

    try {
      const result = await generateText({
        model: this.model,
        ...generateParams
      });

      if (process.env.NODE_ENV === 'development') {
        console.log(`[AIService:${this.config.provider}] @ai-sdk response received`);
      }

      return result;
    } catch (error) {
      console.error(`[AIService:${this.config.provider}] @ai-sdk error:`, error.message);
      throw error;
    }
  }

  /**
   * PARAMETER FORMAT CONVERSION - INTERNAL ↔ SDK
   *
   * Translates between this service's API and @ai-sdk's expectations.
   * Handles message format conversion and parameter mapping.
   *
   * INPUT: Internal format (messages/input, max_tokens, temperature, etc.)
   * OUTPUT: @ai-sdk format (messages/prompt, maxTokens, temperature, etc.)
   *
   * WHY NEEDED:
   * - Different SDKs expect different parameter names
   * - Message format variations
   * - Maintains clean internal API
   */
  _convertToSDKFormat(params) {
    /**
     * PARAMETER MAPPING (SINGLE LOCATION)
     * - This is the only place we translate internal params → SDK params.
     * - If a new param is needed, add it here (don’t spread params at call sites).
     */
    const sdkParams = {};

    if (params.messages && Array.isArray(params.messages)) {
      // Convert to @ai-sdk messages format
      sdkParams.messages = params.messages.map(msg => ({
        role: msg.role || 'user',
        content: msg.content || ''
      }));
    } else {
      // Use prompt format
      sdkParams.prompt = params.input || params.prompt || '';
    }

    // Add optional parameters if provided
    if (params.max_tokens) sdkParams.maxTokens = params.max_tokens;
    if (params.temperature !== undefined) sdkParams.temperature = params.temperature;
    if (params.top_p !== undefined) sdkParams.topP = params.top_p;
    if (params.stream !== undefined) sdkParams.stream = params.stream;
    
    // Agentic Capabilities
    if (params.tools) sdkParams.tools = params.tools;
    if (params.maxSteps) sdkParams.maxSteps = params.maxSteps;

    return sdkParams;
  }

  /**
   * Get mock response for testing
   *
   * CRITICAL: This is ONLY for testing - NO production mock responses
   * FORBIDDEN: Do NOT use this outside of test environments
   */
  _getMockResponse(params) {
    return {
      text: 'This is a mock response for testing purposes.',
      finishReason: 'stop',
      usage: {
        promptTokens: 10,
        completionTokens: 20,
        totalTokens: 30
      }
    };
  }


  /**
   * RESPONSE NORMALIZATION - SDK ↔ INTERNAL
   *
   * Converts @ai-sdk responses to this service's standard format.
   * Handles JSON parsing and metadata extraction.
   *
   * INPUT: @ai-sdk generateText result
   * OUTPUT: Normalized response object with content, usage, metadata
   *
   * WHY IMPORTANT:
   * - Consistent response format for all callers
   * - Automatic JSON parsing for structured outputs
   * - Unified error handling
   */
  _parseSDKResponse(response) {
    /**
     * RESPONSE NORMALIZATION (SINGLE LOCATION)
     * - Callers expect either parsed JSON or a normalized `{content, ...}` object
     * - Do not parse raw provider responses; only handle `@ai-sdk generateText()` output
     */
    try {
      // CRITICAL: Handle ONLY @ai-sdk generateText response format
      // FORBIDDEN: Do NOT parse raw provider API responses

      // Handle Tool Calls (Agentic Mode)
      if (response.toolCalls && response.toolCalls.length > 0) {
        return {
          content: response.text || '', // Can be empty if just calling tools
          toolCalls: response.toolCalls,
          role: 'assistant',
          finish_reason: response.finishReason,
          usage: response.usage,
          model: this._getCurrentProviderConfig().model
        };
      }

      if (!response?.text) {
        throw new Error('Invalid @ai-sdk response: missing text field');
      }

      const content = response.text.trim();
      if (!content) {
        throw new Error('Empty response content from AI');
      }

      // Try to parse as JSON first
      try {
        return JSON.parse(content);
      } catch {
        return {
          content: content,
          role: 'assistant',
          finish_reason: response.finishReason || 'stop',
          usage: response.usage,
          model: this._getCurrentProviderConfig().model
        };
      }
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
   * PROVIDER SWITCHING - RUNTIME RECONFIGURATION
   *
   * Enables changing AI providers without service restart.
   * Reloads all configuration from environment variables.
   *
   * PROCESS:
   * 1. Validate new provider exists
   * 2. Reload all config from .env
   * 3. Reinitialize model with new provider
   *
   * BENEFITS:
   * - Zero-downtime provider switching
   * - Automatic config refresh
   * - Consistent state after switch
   */
  switchProvider(newProvider) {
    if (!this.providerConfigs[newProvider]) {
      throw new Error(`Unsupported provider: ${newProvider}. Supported: ${Object.keys(this.providerConfigs).join(', ')}`);
    }

    // Always reload fresh config from environment
    this.config.provider = newProvider;
    this.providerConfigs = this._createProviderConfigs();

    // Reinitialize with new provider settings
    this.initializeModel();
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
    return ['openai', 'xai', 'google'];
  }

  /**
   * LEGACY PROVIDER DEFAULTS - DO NOT USE FOR AI CONFIGURATION
   *
   * WARNING: This method contains OUTDATED hardcoded values and exists ONLY for legacy compatibility.
   *
   * CRITICAL RESTRICTIONS:
   * - NEVER use for actual LLM configuration
   * - NEVER reference model names from this method
   * - NEVER use for runtime decisions
   *
   * CORRECT APPROACH:
   * - Use _createProviderConfigs() for all AI configuration
   * - Let SDKs handle model versioning automatically
   * - Configure via environment variables only
   *
   * This method will be removed in future versions.
   */
  static getProviderDefaults(provider) {
    // NOTE: These values are intentionally outdated to prevent misuse
    const legacyDefaults = {
      openai: {
        baseURL: 'https://api.openai.com/v1',
        model: 'gpt-3.5-turbo'  // LEGACY - Do not use
      },
      xai: {
        baseURL: 'https://api.x.ai/v1',
        model: 'grok-3'  // LEGACY - Do not use
      }
    };
    return legacyDefaults[provider] || legacyDefaults.openai;
  }
}

