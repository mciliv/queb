const fs = require('fs');
const path = require('path');

class Structuralizer {
  constructor(dependencies = {}) {
    // Required dependencies
    this.aiClient = dependencies.aiClient;
    this.molecularProcessor = dependencies.molecularProcessor;
    this.nameResolver = dependencies.nameResolver;
    this.promptEngine = dependencies.promptEngine;
    this.errorHandler = dependencies.errorHandler;
    this.logger = dependencies.logger;

    // Optional dependencies
    this.cache = dependencies.cache || null;
    
    // Configuration
    this.config = dependencies.config || {
      aiTimeout: 10000,
      maxRetries: 2,
      cacheEnabled: true
    };
    
    // Validate required dependencies
    this._validateDependencies();
  }
  
  /**
   * @private
   */
  // Requires OPENAI_API_KEY to be present in environment and DI to provide aiClient.
  async _detectObjectInImage(imageBase64, x, y) {
    this.logger.info('Detecting object in image', {
      hasCoordinates: x !== undefined && y !== undefined
    });
    
    const prompt = this.promptEngine.generateDetectionPrompt({ x, y });
    
    const response = await this._callAI({
      messages: [{
        role: 'user',
        content: [
          { type: 'text', text: prompt },
          { 
            type: 'image_url', 
            image_url: { 
              url: imageBase64.startsWith('http') 
                ? imageBase64 
                : `data:image/jpeg;base64,${imageBase64}`
            }
          }
        ]
      }],
      max_tokens: 200
    });
    
    const result = this._parseAIResponse(response);
    
    if (!this.promptEngine.validateResponse('detection', result)) {
      throw new Error('Invalid detection response from AI');
    }
    
    return result;
  }
  /**
   * Validate that all required dependencies are provided
   * @private
   */
  _validateDependencies() {
    const required = [
      'molecularProcessor',
      'nameResolver',
      'promptEngine',
      'errorHandler',
      'logger'
    ];
    
    const missing = required.filter(dep => !this[dep]);
    if (missing.length > 0) {
      throw new Error(`Missing required dependencies: ${missing.join(', ')}`);
    }
  }
  

  async chemicals(payload) {
    const startTime = Date.now();
    
    try {
      if (this.cache && this.config.cacheEnabled) {
        const cacheKey = this._getCacheKey(payload);
        const cached = await this.cache.get(cacheKey);
        if (cached) {
          this.logger.info('Cache hit for prediction', { 
            object: payload.object,
            duration: Date.now() - startTime 
          });
          return cached;
        }
      }
      
      const result = await this._performPrediction(payload);
      
      if (this.cache && this.config.cacheEnabled && result) {
        const cacheKey = this._getCacheKey(payload);
        await this.cache.set(cacheKey, result, 300000); // 5 min TTL
      }
      
      this.logger.info('Prediction completed', {
        object: result.object,
        chemicalCount: result.chemicals?.length || 0,
        duration: Date.now() - startTime
      });
      
      return result;
    } catch (error) {
      const handled = await this.errorHandler.handle(error, {
        operation: 'chemicals',
        payload
      });

      this.logger.error('Prediction failed', {
        error: handled,
        duration: Date.now() - startTime
      });

      const err = new Error(handled.message);
      err.code = handled.code;
      err.recoverable = handled.recoverable;
      err.recovery = handled.recovery;
      err.timestamp = handled.timestamp;
      err.context = handled.context;
      throw err;
    }
  }
  
  /**
   * @private
   */
  async _performPrediction(payload) {
    const {
      object: inputObject,
      imageBase64,
      x,
      y,
      lookupMode = 'GPT-5'
    } = payload;
    
    let objectText = inputObject?.trim() || '';
    let recommendedBox = null;
    let reason = null;
    
    // Step 1: Extract object from image if needed
    if (!objectText && imageBase64) {
      if (!this.aiClient) {
        throw new Error('Image prediction requires AI client');
      }
      
      const detection = await this._detectObjectInImage(imageBase64, x, y);
      objectText = detection.object;
      recommendedBox = detection.recommendedBox;
    }
    
    if (!objectText) {
      throw new Error('No object specified for prediction');
    }
    
    // Step 2: Analyze using AI
    let predictionResult;

    switch (lookupMode) {
      case 'ai':
      case 'GPT-5':
        predictionResult = await this._analyzeWithAI(objectText);
        break;

      default:
        throw new Error(`Unknown lookup mode: ${lookupMode}`);
    }
    
    // Step 3: Generate 3D structures
    const molecules = await this._generateStructures(
      predictionResult.chemicals || []
    );
    
    return {
      object: objectText,
      chemicals: molecules,
      recommendedBox,
      reason: predictionResult.reason || reason
    };
  }

  /**
   * @private
   */
  // Core AI analysis of the provided object text.
  // aiClient is injected by DI (see ServiceProvider) and sourced from OPENAI_API_KEY.
  async _analyzeWithAI(objectText) {
    if (!this.aiClient) {
      throw new Error('AI prediction requires AI client');
    }
    
    const prompt = this.promptEngine.generateChemicalPrompt(
      objectText,
      { includeReason: true }
    );
    
    const response = await this._callAI({
      messages: [{ role: 'user', content: prompt }],
      temperature: 1
    });
    
    const result = this._parseAIResponse(response);
    
    if (!this.promptEngine.validateResponse('chemical', result)) {
      throw new Error('Invalid chemical response from AI');
    }
    
    return result;
  }
  
  /**
   * Call AI service with error handling
   * @private
   */
  async _callAI(params) {
    const model = this.config.get ? this.config.get('openai.model') : 'gpt-4o-mini';
    const timeout = this.config.get ? this.config.get('openai.timeout') : 30000;

    return await Promise.race([
      this.aiClient.chat.completions.create({ model, ...params }),
      new Promise((_, reject) =>
        setTimeout(() => reject(new Error('AI request timeout')), timeout)
      )
    ]);
  }
  
  /**
   * Parse AI response
   * @private
   */
  _parseAIResponse(response) {
    let content = '';
    
    try {
      if (!response?.choices?.[0]) {
        throw new Error('Invalid AI response structure');
      }

      const choice = response.choices[0];
      content = choice.message?.content || choice.text || '';

      return JSON.parse(content);
    } catch (error) {
      const repaired = this.promptEngine.repairJSON(content);
      if (repaired) return repaired;

      throw new Error(`Failed to parse AI response: ${error.message}`);
    }
  }
  
  /**
   * Generate 3D structures for molecules
   * @private
   */
  async _generateStructures(chemicals) {
    const results = [];
    
    for (const chemical of chemicals) {
      try {
        let sdfPath = null;
        let status = 'lookup_required';
        
        // Try to generate SDF if we have SMILES
        if (chemical.smiles) {
          sdfPath = await this.molecularProcessor.generateSDF(
            chemical.smiles,
            false // don't overwrite
          );
          status = sdfPath ? 'ok' : 'generation_failed';
        } 
        // Try to resolve by name if no SMILES
        else if (chemical.name) {
          const resolved = await this.nameResolver.resolveName(chemical.name);
          if (resolved.smiles) {
            sdfPath = await this.molecularProcessor.generateSDF(
              resolved.smiles,
              false
            );
            status = sdfPath ? 'ok' : 'generation_failed';
          }
        }
        
        results.push({
          name: chemical.name,
          smiles: chemical.smiles,
          sdfPath,
          status
        });
      } catch (error) {
        this.logger.warn('Failed to generate structure', {
          chemical: chemical.name,
          error: error.message
        });
        
        results.push({
          name: chemical.name,
          smiles: chemical.smiles,
          sdfPath: null,
          status: 'error'
        });
      }
    }
    
    return results;
  }
  
  /**
   * Generate cache key for payload
   * @private
   */
  _getCacheKey(payload) {
    const key = {
      object: payload.object,
      lookupMode: payload.lookupMode,
      x: payload.x,
      y: payload.y,
      imageHash: payload.imageBase64 
        ? this._hashString(payload.imageBase64.substring(0, 100))
        : null
    };
    
    return `structuralizer:${JSON.stringify(key)}`;
  }
  
  /**
   * Simple string hash for cache keys
   * @private
   */
  _hashString(str) {
    let hash = 0;
    for (let i = 0; i < str.length; i++) {
      const char = str.charCodeAt(i);
      hash = ((hash << 5) - hash) + char;
      hash = hash & hash;
    }
    return hash.toString(36);
  }
}

// Export both the class and a legacy wrapper for backwards compatibility
module.exports = Structuralizer;

// Legacy export for existing code
module.exports.chemicals = async function(payload) {
  // This would need to get dependencies from somewhere
  // In practice, this would be removed and all callers would use DI
  throw new Error(
    'Legacy chemicals() function not supported. ' +
    'Please use dependency injection with Structuralizer class.'
  );
};
