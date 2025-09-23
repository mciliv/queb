/**
 * Unified Prompt Engine
 * 
 * This module applies the "Deep Modules" principle by providing a simple interface
 * for prompt generation while hiding the complexity of different prompt types,
 * validation, and context-specific optimizations.
 * 
 * Philosophy: "Information hiding reduces complexity by eliminating dependencies"
 */

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
    // Chemical analysis template
    this._templates.set('chemical_analysis', {
      base: `Analyze the chemical composition of the provided object.
Return JSON only. Provide valid, verified SMILES for each molecule.

RULES:
- Output valid, verified SMILES for each specific molecule
- Prefer canonical, database-standard forms (PubChem/ChEBI)
- Keep SMILES realistic and chemically plausible
- Do not output chemical formulas; always use SMILES
- Avoid commentary; return JSON only
- Focus on characteristic constituents (major first, then minor)
- If nothing identifiable, return empty list but valid JSON

RESPONSE FORMAT:
{
  "object": "Object name",
  "reason": "Brief explanation",
  "chemicals": [
    { "name": "Specific chemical name", "smiles": "valid SMILES" }
  ]
}`,
      examples: this._getChemicalExamples(),
      validation: (result) => this._validateChemicalResponse(result)
    });

    // Object detection template
    this._templates.set('object_detection', {
      base: `Identify the main object in this image and provide a bounding box.
Focus on the most prominent physical object that could contain molecules.

Return JSON format:
{
  "object": "descriptive name",
  "reason": "why this object was identified",
  "recommendedBox": {
    "x": number,
    "y": number, 
    "width": number,
    "height": number
  }
}`,
      validation: (result) => this._validateDetectionResponse(result)
    });

    // Name resolution template
    this._templates.set('name_resolution', {
      base: `Convert object description to specific chemical molecule names.
Output ONLY fully specific molecule names (common names or IUPAC).
No classes/families, no brand names.

Schema:
{
  "object": "string",
  "molecules": [
    {"name": "string", "cid": "number|null"}
  ]
}`,
      validation: (result) => this._validateNameResponse(result)
    });
  }

  /**
   * Get chemical analysis examples based on context
   */
  _getChemicalExamples() {
    return `
EXAMPLES:
// Simple molecules
Water → { "name": "Water", "smiles": "O" }
Ethanol → { "name": "Ethanol", "smiles": "CCO" }
Glucose → { "name": "Glucose", "smiles": "C(C(C(C(C(C=O)O)O)O)O)O" }

// Ionic compounds
Sodium chloride → { "name": "Sodium chloride", "smiles": "[Na+].[Cl-]" }

// Complex materials
Coffee → [
  { "name": "Water", "smiles": "O" },
  { "name": "Caffeine", "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" }
]

Wine → [
  { "name": "Ethanol", "smiles": "CCO" },
  { "name": "Water", "smiles": "O" },
  { "name": "Tartaric acid", "smiles": "OC(C(O)C(O)=O)C(O)=O" }
]`;
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
      (mol.cid === null || typeof mol.cid === 'number')
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

  /**
   * Get context-specific examples
   */
  _getContextExamples(objectType) {
    const contexts = {
      beverage: `
Beverage examples:
- Wine: Ethanol "CCO", Water "O", Tartaric acid
- Coffee: Water "O", Caffeine "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
- Beer: Ethanol "CCO", Water "O", various sugars`,
      
      food: `
Food examples:
- Apple: Water "O", Fructose, Malic acid "C(C(=O)O)C(C(=O)O)O"
- Bread: Starch units, Water "O", various proteins`,
      
      material: `
Material examples:
- Plastic: Polymer units with specific SMILES
- Metal: Ionic forms like "[Na+]" for sodium
- Glass: Silicon dioxide components`
    };

    return contexts[objectType] || '';
  }

  /**
   * PUBLIC INTERFACE - Simple methods hiding complexity
   */

  /**
   * Generate a prompt for chemical analysis
   */
  generateChemicalPrompt(objectDescription, options = {}) {
    const template = this._templates.get('chemical_analysis');
    const contextExamples = this._getContextExamples(options.objectType || 'general');
    
    return [
      template.base,
      `Object: ${objectDescription}`,
      template.examples,
      contextExamples,
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


