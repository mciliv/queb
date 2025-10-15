const fs = require('fs');
const path = require('path');

class PromptEngine {
  constructor() {
    this._templates = new Map();
    this._contexts = new Map();
    this._initializeTemplates();
  }

  /**
   * Initialize prompt templates - complex logic hidden from users
   */
  _initializeTemplates() {
    // Load prompt templates from plain text files (synchronous for simplicity)
    const chemicalBase = this._loadPromptFile('chemical_analysis.txt');
    const detectionBase = this._loadPromptFile('object_detection.txt');
    const nameBase = this._loadPromptFile('name_resolution.txt');

    this._templates.set('chemical_analysis', {
      base: chemicalBase,
      validation: (result) => this._validateChemicalResponse(result)
    });

    this._templates.set('object_detection', {
      base: detectionBase,
      validation: (result) => this._validateDetectionResponse(result)
    });

    this._templates.set('name_resolution', {
      base: nameBase,
      validation: (result) => this._validateNameResponse(result)
    });
  }

<<<<<<< Updated upstream
  /**
   * Get chemical analysis examples based on context
   */
  _getChemicalExamples() {
    return `
GUIDANCE:
- Use authoritative databases: PubChem (NCBI) and ChEBI.
- Prefer PubChem IsomericSMILES when a CID is available.
- Do not guess: if uncertain or ambiguous, set "smiles": null.
- Output JSON only; no commentary.
- Focus on characteristic constituents (major first, then minor).`;
=======
  _loadPromptFile(filename) {
    const filePath = path.join(__dirname, 'prompts', filename);
    try {
      const text = fs.readFileSync(filePath, 'utf8');
      return typeof text === 'string' ? text.trim() : '';
    } catch (err) {
      // If missing, return empty string to avoid runtime failure
      return '';
    }
>>>>>>> Stashed changes
  }

  /**
   * Response validation methods - encapsulate complexity
   */
  _validateChemicalResponse(response) {
    if (!response || typeof response !== 'object') return false;
    if (!response.object || typeof response.object !== 'string') return false;
    if (!Array.isArray(response.chemicals)) return false;
    
    return response.chemicals.every(chem => 
      chem.name && typeof chem.name === 'string' &&
      chem.smiles && typeof chem.smiles === 'string' &&
      this._isValidSMILES(chem.smiles)
    );
  }

  _validateDetectionResponse(response) {
    if (!response || typeof response !== 'object') return false;
    if (!response.object || typeof response.object !== 'string') return false;
    
    // Bounding box is optional but if present must be valid
    if (response.recommendedBox) {
      const box = response.recommendedBox;
      return typeof box.x === 'number' && typeof box.y === 'number' &&
             typeof box.width === 'number' && typeof box.height === 'number' &&
             box.width > 0 && box.height > 0;
    }
    
    return true;
  }

  _validateNameResponse(response) {
    if (!response || typeof response !== 'object') return false;
    if (!response.object || typeof response.object !== 'string') return false;
    if (!Array.isArray(response.molecules)) return false;
    
    return response.molecules.every(mol => 
      mol.name && typeof mol.name === 'string' &&
      true
    );
  }

  /**
   * Basic SMILES validation
   */
  _isValidSMILES(smiles) {
    if (!smiles || typeof smiles !== 'string') return false;
    const cleaned = smiles.trim();
    if (cleaned === '' || cleaned === 'N/A') return false;
    
    // Basic pattern validation
    return /^[A-Za-z0-9\[\]()=+\-#@\\/.]+$/.test(cleaned);
  }

<<<<<<< Updated upstream
  /**
   * Get context-specific examples
   */
  _getContextExamples(objectType) {
    const contexts = {
      beverage: `
Beverage guidance:
- Identify typical constituents for the beverage.
- Verify identities against PubChem/ChEBI; prefer PubChem IsomericSMILES.`,
      
      food: `
Food guidance:
- Identify characteristic small molecules (sugars, acids, etc.).
- Verify against PubChem/ChEBI; return SMILES or null if uncertain.`,
      
      material: `
Material guidance:
- Identify representative units/ions for the material.
- Verify against PubChem/ChEBI; prefer canonical SMILES.`
    };

    return contexts[objectType] || '';
  }
=======
  // _getContextExamples removed: prompts now live in plain text files
>>>>>>> Stashed changes

  /**
   * PUBLIC INTERFACE - Simple methods hiding complexity
   */

  /**
   * Generate a prompt for chemical analysis
   */
  generateChemicalPrompt(objectDescription, options = {}) {
    const template = this._templates.get('chemical_analysis');
    
    // Minimal prompt: keep rules and schema only, avoid extra examples/context that can distract the model.
    // Preserve option to include a brief reason.
    return [
      template.base,
      `Object: ${objectDescription}`,
      options.includeReason ? 'Include a brief reason for your identification.' : ''
    ].filter(Boolean).join('\n\n');
  }

  /**
   * Generate a prompt for object detection from images
   */
  generateDetectionPrompt(coordinates = null) {
    const template = this._templates.get('object_detection');
    let prompt = template.base;
    
    if (coordinates && coordinates.x !== null && coordinates.y !== null) {
      prompt += `\n\nUser clicked at coordinates (${coordinates.x}, ${coordinates.y}). Focus analysis around this area.`;
    }
    
    return prompt;
  }

  /**
   * Generate a prompt for name resolution
   */
  generateNamePrompt(objectDescription) {
    const template = this._templates.get('name_resolution');
    
    return [
      template.base,
      `Object: ${objectDescription}`,
      `Order molecules by abundance in the object.`,
      `Prefer specific isomers when typical (e.g., "L-ascorbic acid", "(+)-catechin").`
    ].join('\n\n');
  }

  /**
   * Validate AI response against expected format
   */
  validateResponse(promptType, response) {
    const template = this._templates.get(promptType);
    if (!template || !template.validation) return true;
    
    try {
      const parsed = typeof response === 'string' ? JSON.parse(response) : response;
      return template.validation(parsed);
    } catch (error) {
      return false;
    }
  }

  /**
   * Get available prompt types
   */
  getAvailablePrompts() {
    return Array.from(this._templates.keys());
  }

  /**
   * Repair incomplete JSON responses
   */
  repairJSON(jsonString) {
    if (!jsonString || typeof jsonString !== 'string') return null;
    
    try {
      return JSON.parse(jsonString);
    } catch (error) {
      // Try to extract and repair JSON
      const match = jsonString.match(/\{[\s\S]*\}/);
      if (!match) return null;
      
      let jsonStr = match[0].trim();
      
      // Count and balance braces/brackets
      const openBraces = (jsonStr.match(/\{/g) || []).length;
      const closeBraces = (jsonStr.match(/\}/g) || []).length;
      const openBrackets = (jsonStr.match(/\[/g) || []).length;
      const closeBrackets = (jsonStr.match(/\]/g) || []).length;
      
      // Add missing closing brackets/braces
      for (let i = 0; i < openBrackets - closeBrackets; i++) {
        jsonStr += ']';
      }
      for (let i = 0; i < openBraces - closeBraces; i++) {
        jsonStr += '}';
      }
      
      try {
        return JSON.parse(jsonStr);
      } catch (repairError) {
        return null;
      }
    }
  }
}

// Export singleton instance
const promptEngine = new PromptEngine();

module.exports = promptEngine;


