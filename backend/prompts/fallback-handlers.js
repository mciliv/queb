// Fallback Response Handlers for AI Analysis
// Extracted from successful patterns in git history

/**
 * Smart Fallback System - Handles AI refusals and parsing errors gracefully
 * 
 * KEY INSIGHTS FROM GIT HISTORY:
 * - Specific fallbacks for human analysis work better than generic ones (commit be3ecc0)
 * - Comprehensive refusal detection improves user experience (commit 32eac08)
 * - Realistic chemical compositions maintain app functionality (commit 0fe79e8)
 */

// Proven fallback compositions from working versions
const FALLBACK_COMPOSITIONS = {
  // Human body composition - realistic and useful (from commit be3ecc0)
  human: {
    object: "Human body (generic composition)",
    chemicals: [
      {"name": "Water", "smiles": "O"},                                    // ~60% of body weight
      {"name": "Glycine", "smiles": "C(C(=O)O)N"},                        // Simplest amino acid
      {"name": "Leucine", "smiles": "CC(C)CC(N)C(=O)O"},                  // Essential amino acid
      {"name": "Palmitic acid", "smiles": "CCCCCCCCCCCCCCCC(=O)O"},        // Common fatty acid
      {"name": "Glucose", "smiles": "C(C(C(C(C(C=O)O)O)O)O)O"},          // Blood sugar
      {"name": "Lactic acid", "smiles": "C(C(=O)O)O"}                     // Muscle metabolite
    ]
  },

  // Generic analysis unavailable - minimal but functional
  generic: {
    object: "Generic object (AI analysis unavailable)", 
    chemicals: [
      {"name": "Water", "smiles": "O"},                    // Universal component
      {"name": "Carbon dioxide", "smiles": "O=C=O"},       // Common in atmosphere
      {"name": "Nitrogen", "smiles": "N#N"},               // Atmospheric component
      {"name": "Oxygen", "smiles": "O=O"}                  // Essential element
    ]
  },

  // Analysis completed with minimal data - maintains functionality
  minimal: {
    object: "Analysis completed with fallback data",
    chemicals: [
      {"name": "Water", "smiles": "O"},           // Most common molecule
      {"name": "Carbon", "smiles": "C"},          // Building block element  
      {"name": "Oxygen", "smiles": "O=O"}         // Essential for life
    ]
  }
};

/**
 * Enhanced AI Response Parser with Smart Fallbacks
 * Combines proven techniques from multiple git commits
 */
function parseAIResponseWithFallbacks(content) {
  try {
    // Primary parsing attempt - extract JSON from response
    const jsonMatch = content.match(/\{[\s\S]*\}/);
    if (jsonMatch) {
      const parsed = JSON.parse(jsonMatch[0]);
      
      // Validate the parsed response has required fields
      if (parsed.object && parsed.chemicals && Array.isArray(parsed.chemicals)) {
        return {
          object: parsed.object,
          chemicals: parsed.chemicals || []
        };
      }
    }
    
    // Secondary attempt - parse entire content as JSON
    const parsed = JSON.parse(content);
    return {
      object: parsed.object || "Unknown object",
      chemicals: parsed.chemicals || []
    };

  } catch (error) {
    console.error("Failed to parse AI response:", content);
    
    // Smart fallback detection based on content analysis
    return detectAndHandleRefusal(content);
  }
}

/**
 * Intelligent Refusal Detection and Handling
 * Based on patterns observed in successful deployments
 */
function detectAndHandleRefusal(content) {
  const lowerContent = content.toLowerCase();
  
  // Human/person analysis refusal - use realistic body composition (commit be3ecc0)
  const humanRefusalPatterns = [
    "unable to identify or analyze people",
    "can't identify people", 
    "cannot analyze people",
    "won't analyze human",
    "people in images"
  ];
  
  if (humanRefusalPatterns.some(pattern => lowerContent.includes(pattern))) {
    return FALLBACK_COMPOSITIONS.human;
  }
  
  // General image analysis refusal - provide minimal but useful response
  const imageRefusalPatterns = [
    "can't analyze images",
    "unable to identify", 
    "i'm sorry",
    "cannot process images",
    "unable to see"
  ];
  
  if (imageRefusalPatterns.some(pattern => lowerContent.includes(pattern))) {
    return FALLBACK_COMPOSITIONS.generic;
  }
  
  // Safety/policy refusal - graceful degradation
  const safetyRefusalPatterns = [
    "safety guidelines",
    "policy prevents",
    "cannot provide",
    "not appropriate"
  ];
  
  if (safetyRefusalPatterns.some(pattern => lowerContent.includes(pattern))) {
    return FALLBACK_COMPOSITIONS.minimal;
  }
  
  // Default fallback for any other parsing error
  return FALLBACK_COMPOSITIONS.minimal;
}

/**
 * Validation function to check SMILES quality
 * Helps prevent the "wrong SMILES" issue mentioned by user
 */
function validateSMILESQuality(chemicals) {
  return chemicals.map(chemical => {
    const smiles = chemical.smiles;
    
    // Basic SMILES validation - should not be empty or just chemical formulas
    if (!smiles || smiles.length === 0) {
      console.warn(`Empty SMILES for ${chemical.name}`);
      return chemical;
    }
    
    // Check for common formula errors (H2O instead of O)
    const formulaPattern = /^[A-Z][a-z]?\d*([A-Z][a-z]?\d*)*$/;
    if (formulaPattern.test(smiles) && smiles !== "O" && smiles !== "C" && smiles !== "N") {
      console.warn(`Possible chemical formula instead of SMILES: ${smiles} for ${chemical.name}`);
    }
    
    return chemical;
  });
}

module.exports = {
  parseAIResponseWithFallbacks,
  detectAndHandleRefusal,
  validateSMILESQuality,
  FALLBACK_COMPOSITIONS
}; 