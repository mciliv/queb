/**
 * Molecular Prediction Service
 * 
 * This is a deep module that provides a simple interface for molecular prediction
 * while hiding the complexity of AI integration, data processing, and error handling.
 * 
 * Philosophy: "The most important issue in software design is complexity" - John Ousterhout
 * This service reduces complexity by providing a unified interface for all molecular prediction tasks.
 */

const { config } = require('./services');
const promptEngine = require('./PromptEngine');
const fs = require('fs');
const path = require('path');

class MolecularPredictionService {
  constructor(dependencies = {}) {
    // Dependency injection for testability
    this.aiClient = dependencies.aiClient;
    this.molecularProcessor = dependencies.molecularProcessor;
    this.nameResolver = dependencies.nameResolver;
    this.logger = dependencies.logger || console;
    
    // Internal state
    this._isInitialized = false;
    this._capabilities = new Set();
  }

  /**
   * Initialize the service with required dependencies
   * Deep module: complex initialization logic hidden from users
   */
  async initialize() {
    if (this._isInitialized) return;

    try {
      // Load dependencies dynamically to avoid circular imports
      if (!this.aiClient) {
        const { createOpenAIClient } = require('./services');
        this.aiClient = createOpenAIClient();
      }

      if (!this.molecularProcessor) {
        const MolecularProcessor = require('../server/services/molecular-processor');
        this.molecularProcessor = new MolecularProcessor();
      }

      if (!this.nameResolver) {
        this.nameResolver = require('../server/services/name-resolver');
      }

      // Determine capabilities based on configuration
      const openaiKey = config.get('openai.apiKey');
      if (this.aiClient && openaiKey) {
        this._capabilities.add('ai_prediction');
        this._capabilities.add('image_prediction');
      }

      this._capabilities.add('molecular_processing');
      this._capabilities.add('name_resolution');

      this._isInitialized = true;
      this.logger.info('✅ Molecular Prediction Service initialized');
      
    } catch (error) {
      this.logger.error('❌ Failed to initialize Molecular Prediction Service:', error);
      throw new Error('Service initialization failed');
    }
  }

  /**
   * Predict molecular information from text description
   * Simple interface hiding complex AI interaction and data processing
   */
  async predictText(description, options = {}) {
    await this.initialize();
    
    if (!description || typeof description !== 'string') {
      throw new Error('Invalid description: must be a non-empty string');
    }

    if (!this._capabilities.has('ai_prediction')) {
      throw new Error('AI prediction not available - check API key configuration');
    }

    try {
      this.logger.info(`🔍 Predicting from text: "${description}"`);

      // Load minimal chemical prediction prompt text + optional reason suffix and compose inline
      let chemBase = '';
      let reasonSuffix = '';
      try {
        const chemPath = path.join(__dirname, 'prompts', 'chemical_analysis.txt');
        chemBase = fs.readFileSync(chemPath, 'utf8').trim();
      } catch (_) { chemBase = ''; }
      try {
        const reasonPath = path.join(__dirname, 'prompts', 'reason_suffix.txt');
        reasonSuffix = fs.readFileSync(reasonPath, 'utf8').trim();
      } catch (_) { reasonSuffix = ''; }
      const prompt = [
        chemBase,
        `Object: ${description}`,
        (options.includeReason !== false) ? reasonSuffix : ''
      ].filter(Boolean).join('\n\n');

      // Call AI service with error handling
      const response = await this._callAI({
        messages: [{ role: 'user', content: prompt }],
        response_format: { type: 'json_object' }
      });

      // Parse and validate response
      const result = this._parseAIResponse(response);
      
      // Enrich with additional data
      const enriched = await this._enrichMolecularData(result.chemicals || []);
      
      // Generate 3D structures
      const withStructures = await this._generateStructures(enriched);

      this.logger.info(`✅ Prediction complete: Found ${withStructures.length} molecules`);

      return {
        object: result.object || description,
        reason: result.reason,
        molecules: withStructures,
        metadata: {
          processingTime: Date.now(),
          capabilities: Array.from(this._capabilities),
          enrichmentApplied: true
        }
      };

    } catch (error) {
      this.logger.error('❌ Text prediction failed:', error);
      throw this._createPredictionError(error, 'text_prediction');
    }
  }

  /**
   * Predict molecular information from image
   * Simple interface hiding complex image processing and AI vision
   */
  async predictImage(imageData, coordinates = null, options = {}) {
    await this.initialize();

    if (!this._capabilities.has('image_prediction')) {
      throw new Error('Image prediction not available - check AI service configuration');
    }

    try {
      this.logger.info('🖼️ Predicting from image for molecular content');

      // First, detect the main object in the image
      const detection = await this._detectObject(imageData, coordinates);
      
      // Then predict the detected object for molecular content
      if (detection.object && detection.object !== 'Unknown object') {
        const textResult = await this.predictText(detection.object, {
          ...options,
          objectType: this._detectObjectType(detection.object)
        });

        return {
          ...textResult,
          imagePrediction: {
            detectedObject: detection.object,
            reason: detection.reason,
            recommendedBox: detection.recommendedBox
          }
        };
      }

      // Fallback for unknown objects
      return {
        object: 'Unknown object',
        reason: 'Could not identify a specific object for molecular prediction',
        molecules: [],
        imagePrediction: detection
      };

    } catch (error) {
      this.logger.error('❌ Image prediction failed:', error);
      throw this._createPredictionError(error, 'image_prediction');
    }
  }

  /**
   * Generate 3D molecular structures from SMILES
   * Simple interface hiding complex molecular processing
   */
  async generateStructures(smilesArray, options = {}) {
    await this.initialize();

    if (!Array.isArray(smilesArray)) {
      throw new Error('Invalid input: smilesArray must be an array');
    }

    try {
      const result = await this.molecularProcessor.processSmiles(smilesArray, options.overwrite);
      
      return {
        structures: result.sdfPaths,
        errors: result.errors,
        skipped: result.skipped,
        summary: `Generated ${result.sdfPaths.length} structures from ${smilesArray.length} SMILES`
      };

    } catch (error) {
      this.logger.error('❌ Structure generation failed:', error);
      throw this._createPredictionError(error, 'structure_generation');
    }
  }

  /**
   * PRIVATE METHODS - Implementation details hidden from users
   */

  /**
   * Call AI service with standardized error handling
   */
  async _callAI(params) {
    const model = config.get('openai.model');
    const timeout = config.get('openai.timeout');

    return await Promise.race([
      this.aiClient.chat.completions.create({ model, ...params }),
      new Promise((_, reject) => 
        setTimeout(() => reject(new Error('AI request timeout')), timeout)
      )
    ]);
  }

  /**
   * Parse and validate AI response
   */
  _parseAIResponse(response) {
    const content = response?.choices?.[0]?.message?.content;
    if (!content) throw new Error('Empty AI response');

    // Try parsing with repair if needed
    let parsed = promptEngine.repairJSON(content);
    if (!parsed) {
      throw new Error('Could not parse AI response as JSON');
    }

    return parsed;
  }

  /**
   * Detect object type for context-aware processing
   */
  _detectObjectType(description) {
    const lower = description.toLowerCase();
    
    if (lower.includes('wine') || lower.includes('beer') || lower.includes('coffee') || 
        lower.includes('drink') || lower.includes('beverage')) {
      return 'beverage';
    }
    
    if (lower.includes('fruit') || lower.includes('food') || lower.includes('apple') || 
        lower.includes('vegetable')) {
      return 'food';
    }
    
    if (lower.includes('plastic') || lower.includes('metal') || lower.includes('wood')) {
      return 'material';
    }
    
    return 'general';
  }

  /**
   * Detect object in image using AI vision
   */
  async _detectObject(imageData, coordinates) {
    // Load detection prompt base text and append coordinates if provided
    let detectionBase = '';
    try {
      const detPath = path.join(__dirname, 'prompts', 'object_detection.txt');
      detectionBase = fs.readFileSync(detPath, 'utf8').trim();
    } catch (_) { detectionBase = ''; }
    let prompt = detectionBase;
    if (coordinates && Number.isFinite(coordinates.x) && Number.isFinite(coordinates.y)) {
      prompt += `\n\nUser clicked at coordinates (${coordinates.x}, ${coordinates.y}). Focus prediction around this area.`;
    }

    const messages = [{
      role: 'user',
      content: [
        { type: 'text', text: prompt },
        { type: 'image_url', image_url: { url: imageData, detail: 'high' } }
      ]
    }];

    const response = await this._callAI({
      messages,
      max_tokens: 500,
      response_format: { type: 'json_object' }
    });

    return this._parseAIResponse(response);
  }

  /**
   * Enrich molecular data with additional information
   */
  async _enrichMolecularData(molecules) {
    const enriched = [];
    
    for (const molecule of molecules) {
      let enrichedMol = { ...molecule };
      
      // Try to resolve missing SMILES or additional properties
      if (!enrichedMol.smiles && enrichedMol.name) {
        try {
          const resolved = await this.nameResolver.resolveName(enrichedMol.name);
          if (resolved) {
            enrichedMol.smiles = resolved.smiles || enrichedMol.smiles;
            enrichedMol.name = resolved.title || resolved.iupac || enrichedMol.name;
          }
        } catch (error) {
          // Continue without enrichment
          this.logger.warn(`⚠️ Could not enrich molecule: ${enrichedMol.name}`);
        }
      }
      
      enriched.push(enrichedMol);
    }
    
    return enriched;
  }

  /**
   * Generate 3D structures for molecules
   */
  async _generateStructures(molecules) {
    const withStructures = [];
    
    for (const molecule of molecules) {
      let structurePath = null;
      
      try {
        if (molecule.smiles) {
          structurePath = await this.molecularProcessor.generateSDF(molecule.smiles, false);
        } else if (molecule.name) {
          const byName = await this.molecularProcessor.generateSDFByName(molecule.name, false);
          structurePath = byName?.sdfPath || null;
        }
      } catch (error) {
        this.logger.warn(`⚠️ Structure generation failed for ${molecule.name}:`, error.message);
      }
      
      withStructures.push({
        ...molecule,
        sdfPath: structurePath,
        status: structurePath ? 'ok' : (molecule.smiles ? 'smiles_only' : 'lookup_required')
      });
    }
    
    return withStructures;
  }

  /**
   * Create standardized error objects
   */
  _createPredictionError(originalError, context) {
    const error = new Error(`Molecular prediction failed: ${originalError.message}`);
    error.context = context;
    error.originalError = originalError;
    error.timestamp = new Date().toISOString();
    return error;
  }

  /**
   * Get service status and capabilities
   */
  getStatus() {
    return {
      initialized: this._isInitialized,
      capabilities: Array.from(this._capabilities),
      configuration: config.getDebugInfo()
    };
  }
}

// Export singleton instance
const molecularPredictionService = new MolecularPredictionService();

module.exports = molecularPredictionService;


