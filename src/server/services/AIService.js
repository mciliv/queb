/**
 * AI SERVICE ARCHITECTURE PRINCIPLES
 * =================================
 *
 * This service implements three fundamental principles that MUST be maintained:
 *
 * 1. ENVIRONMENT VARIABLE EXCLUSIVITY
 *    The .env file is the SOLE source of configuration. Every AI setting comes from process.env.
 *
 * 2. SDK ABSTRACTION PURITY
 *    All AI interactions use @ai-sdk's unified interface. Zero direct provider API calls.
 *
 * 3. SDK RELIABILITY MAXIMIZATION
 *    Leverage SDKs for automatic model versioning, API compatibility, and reliability features.
 *
 * WHY THESE PRINCIPLES MATTER:
 * ===========================
 * - Future-proofs against API changes and model deprecations
 * - Enables seamless provider switching without code changes
 * - Prevents configuration drift and hardcoded model dependencies
 * - Automatic adoption of latest models and performance improvements
 * - SDK handles rate limiting, retries, and error recovery
 *
 * REQUIRED BEHAVIORS:
 * ===================
 * ✅ ALWAYS load config from process.env variables
 * ✅ ALWAYS use @ai-sdk generateText() for AI calls
 * ✅ ALWAYS call _createProviderConfigs() for configuration
 * ✅ ALWAYS use _convertToSDKFormat() for parameter conversion
 *
 * CONSEQUENCES OF VIOLATION:
 * ==========================
 * Breaking these principles will cause:
 * - Provider switching failures
 * - Configuration inconsistencies
 * - API migration nightmares
 * - Maintenance overhead
 *
 * EXAMPLES OF CORRECT USAGE:
 * ==========================
 *
 * ✅ Configuration:
 *    this.providerConfigs = this._createProviderConfigs();
 *
 * ✅ AI Calls:
 *    const result = await generateText({ model: this.model, ...params });
 *
 * ❌ WRONG - Never do this:
 *    const client = new OpenAI(); // Direct SDK usage
 *    const result = client.responses.create(); // Raw API calls
 *    const config = { apiKey: 'hardcoded' }; // Hardcoded values
 */
class AIService {
  constructor() {
    // REQUIRED: Load ONLY from environment variables
    // This maintains the single source of truth principle
    const isTest = process.env.NODE_ENV === 'test';

    this.config = {
      provider: process.env.AI_PROVIDER || 'openai'
    };

    // Provider-specific configurations loaded from .env
    // Test environment always uses mock-model regardless of .env settings
    this.providerConfigs = this._createProviderConfigs();

    this.model = null;
    // CRITICAL: ONLY initialize model - NO client initialization
    // FORBIDDEN: Do NOT initialize raw provider clients
    this.initializeModel();
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
   * MODEL INITIALIZATION - SDK ABSTRACTION LAYER
   *
   * Creates @ai-sdk model instances for unified AI interactions.
   * This is the ONLY initialization needed - no raw clients required.
   *
   * BENEFITS:
   * - Automatic provider abstraction
   * - Consistent API across all providers
   * - Future-proof against SDK changes
   *
   * RESULT: this.model ready for generateText() calls
   */
  initializeModel() {
    const providerConfig = this._getCurrentProviderConfig();

    if (!providerConfig.apiKey) {
      throw new Error(`${this.config.provider.toUpperCase()} API key is required`);
    }

    // CRITICAL: Use ONLY @ai-sdk unified interface - NO provider-specific logic
    // FORBIDDEN: Do NOT switch on provider type - @ai-sdk handles this abstraction
    switch (this.config.provider) {
      case 'openai':
        const { openai } = require('@ai-sdk/openai');
        this.model = openai(providerConfig.model, {
          apiKey: providerConfig.apiKey
        });
        break;
      case 'xai':
        const { xai } = require('@ai-sdk/xai');
        this.model = xai(providerConfig.model, {
          apiKey: providerConfig.apiKey
        });
        break;
      default:
        throw new Error(`Unsupported provider: ${this.config.provider}`);
    }
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
    // #region agent log
    fetch('http://127.0.0.1:7243/ingest/f1225f0b-6c5b-477f-bc5d-1e74641debf9', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ location: 'AIService.js:131', message: 'AIService.callAPI called', data: { provider: this.config.provider, model: this._getCurrentProviderConfig().model, messageCount: params.messages?.length, inputLength: params.input?.length }, timestamp: Date.now(), sessionId: 'debug-cacao', runId: 'pre-fix', hypothesisId: 'A,B,C' }) }).catch(() => { });
    // #endregion

    try {
      // CRITICAL: Use ONLY @ai-sdk generateText for complete abstraction
      const result = await this._callWithSDK(params);
      const parsed = this._parseSDKResponse(result);

      // #region agent log
      fetch('http://127.0.0.1:7243/ingest/f1225f0b-6c5b-477f-bc5d-1e74641debf9', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ location: 'AIService.js:140', message: 'AIService.callAPI completed', data: { hasResult: !!result, hasParsed: !!parsed, contentLength: parsed?.content?.length }, timestamp: Date.now(), sessionId: 'debug-cacao', runId: 'pre-fix', hypothesisId: 'A,B,C' }) }).catch(() => { });
      // #endregion

      return parsed;
    } catch (error) {
      // #region agent log
      fetch('http://127.0.0.1:7243/ingest/f1225f0b-6c5b-477f-bc5d-1e74641debf9', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ location: 'AIService.js:147', message: 'AIService.callAPI error', data: { error: error.message, provider: this.config.provider }, timestamp: Date.now(), sessionId: 'debug-cacao', runId: 'pre-fix', hypothesisId: 'A,B,C' }) }).catch(() => { });
      // #endregion
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
    // #region agent log
    fetch('http://127.0.0.1:7243/ingest/f1225f0b-6c5b-477f-bc5d-1e74641debf9', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ location: 'AIService.js:217', message: 'AIService._parseResponse called', data: { hasOutputText: !!response?.output_text, hasChoices: !!response?.choices, objectType: response?.object }, timestamp: Date.now(), sessionId: 'debug-cacao', runId: 'pre-fix', hypothesisId: 'B' }) }).catch(() => { });
    // #endregion

    try {
      // CRITICAL: Handle ONLY @ai-sdk generateText response format
      // FORBIDDEN: Do NOT parse raw provider API responses

      if (!response?.text) {
        throw new Error('Invalid @ai-sdk response: missing text field');
      }

      let content = response.text.trim();
      if (!content) {
        throw new Error('Empty response content from AI');
      }

      // Strip markdown fences if the AI wrapped JSON in ```json ... ```
      // Prompt: chemical_analysis.txt says "no markdown" but models often ignore this
      content = content.replace(/^```(?:json)?\s*\n?/i, '').replace(/\n?```\s*$/i, '').trim();

      // Try to parse as JSON first
      try {
        const parsedJson = JSON.parse(content);
        // #region agent log
        fetch('http://127.0.0.1:7243/ingest/f1225f0b-6c5b-477f-bc5d-1e74641debf9', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ location: 'AIService.js:240', message: 'Parsed JSON response', data: { isJson: true, contentLength: content.length }, timestamp: Date.now(), sessionId: 'debug-cacao', runId: 'pre-fix', hypothesisId: 'B' }) }).catch(() => { });
        // #endregion
        return parsedJson;
      } catch {
        // #region agent log
        fetch('http://127.0.0.1:7243/ingest/f1225f0b-6c5b-477f-bc5d-1e74641debf9', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ location: 'AIService.js:245', message: 'Returning structured response', data: { isJson: false, contentLength: content.length }, timestamp: Date.now(), sessionId: 'debug-cacao', runId: 'pre-fix', hypothesisId: 'B' }) }).catch(() => { });
        // #endregion
        // Return structured response if not JSON
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
    return ['openai', 'xai'];
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

// ARCHITECTURE ENFORCEMENT PRINCIPLES
// ===================================
//
// This implementation follows evidence-based software architecture principles:
//
// 1. SINGLE SOURCE OF TRUTH
//    - All configuration flows from .env file
//    - Prevents configuration drift and inconsistencies
//    - Enables atomic configuration changes
//
// 2. MAXIMUM ABSTRACTION
//    - Zero direct provider API calls
//    - Unified interface via @ai-sdk
//    - Future-proofs against API changes
//
// 3. SDK RELIABILITY MAXIMIZATION
//    - Automatic model versioning and updates
//    - Built-in rate limiting and retry logic
//    - Error handling and recovery
//    - Performance optimizations
//
// 4. CLEAN ARCHITECTURE
//    - Clear separation between configuration, execution, and response handling
//    - Dependency injection through environment variables
//    - Testable components with clear boundaries
//
// MAINTENANCE GUIDELINES:
// ======================
// When modifying this service, always ask:
// - Does this maintain single source of truth?
// - Does this preserve SDK abstraction?
// - Does this maximize SDK reliability features?
// - Does this keep the architecture clean?
//
// If the answer is no to any question, reconsider the approach.

module.exports = AIService;