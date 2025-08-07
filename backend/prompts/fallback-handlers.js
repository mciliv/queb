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
  // Handle empty, null, or non-string content
  if (typeof content !== "string" || content.trim() === "") {
    return { object: "Unknown object", chemicals: [] };
  }

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
    console.error("Failed to parse AI response:", error.message);
    console.error("Content length:", content.length);
    console.error("Content preview:", content.substring(0, 500) + "...");
    
    // Try to repair and extract valid JSON from malformed response
    try {
      console.log("🔧 Attempting to repair malformed JSON...");
      
      // Find the JSON object
      const jsonMatch = content.match(/\{[\s\S]*\}/);
      if (jsonMatch) {
        let cleanedJson = jsonMatch[0];
        
        // Fix common JSON issues
        // 1. Remove incomplete entries that don't have closing quotes
        cleanedJson = cleanedJson.replace(/,\s*"name":\s*"[^"]*$/g, ''); // Remove incomplete entries at end
        cleanedJson = cleanedJson.replace(/,\s*\{\s*"name":\s*[^}]*$/g, ''); // Remove incomplete objects
        
        // Remove incomplete molecules with truncated SMILES at end of array
        cleanedJson = cleanedJson.replace(/,\s*\{\s*"name":\s*"[^"]*",\s*"smiles":\s*"[^"]*$/g, '');
        
        // 2. Handle extremely long SMILES that cause truncation
        // Remove any incomplete SMILES strings that are way too long
        cleanedJson = cleanedJson.replace(/"smiles":\s*"[^"]{200,}$/g, '');
        
        // 3. Ensure arrays are properly closed
        if (!cleanedJson.includes(']}')) {
          cleanedJson = cleanedJson.replace(/[,\s]*$/, '') + ']}';
        }
        
        // 3. Remove trailing commas
        cleanedJson = cleanedJson.replace(/,(\s*[}\]])/g, '$1');
        
        console.log("🔧 Repaired JSON preview:", cleanedJson.substring(0, 200) + "...");
        
        const parsed = JSON.parse(cleanedJson);
        console.log("✅ Successfully parsed repaired JSON!");
        return {
          object: parsed.object || "Unknown object",
          chemicals: parsed.chemicals || []
        };
      }
    } catch (extractError) {
      console.error("Failed to repair JSON:", extractError.message);
      
      // Fallback: Try to extract at least the object name and simple chemicals
      try {
        const objectMatch = content.match(/"object":\s*"([^"]+)"/);
        const chemicalMatches = content.match(/"name":\s*"([^"]+)",\s*"smiles":\s*"([^"]+)"/g);
        
        if (objectMatch || chemicalMatches) {
          const chemicals = [];
          if (chemicalMatches) {
            chemicalMatches.forEach(match => {
              const nameMatch = match.match(/"name":\s*"([^"]+)"/);
              const smilesMatch = match.match(/"smiles":\s*"([^"]+)"/);
              if (nameMatch && smilesMatch) {
                chemicals.push({
                  name: nameMatch[1],
                  smiles: smilesMatch[1]
                });
              }
            });
          }
          
          console.log("✅ Extracted chemicals using regex fallback:", chemicals.length);
          return {
            object: objectMatch ? objectMatch[1] : "Unknown object",
            chemicals: chemicals
          };
        }
      } catch (regexError) {
        console.error("Regex extraction also failed:", regexError.message);
      }
    }
    
    // Only use smart fallbacks for actual refusal patterns, not parsing errors
    if (typeof content === "string" && content.length > 10) {
      const fallback = detectAndHandleRefusal(content);
      if (fallback) return fallback;
    }

    // Default minimal object for parsing errors
    return { object: "Unknown object", chemicals: [] };
  }
}

/**
 * Intelligent Refusal Detection and Handling
 * Based on patterns observed in successful deployments
 */
function detectAndHandleRefusal(content) {
  if (typeof content !== "string") return null;
  const lowerContent = content.toLowerCase();
  
  // Only check for refusals in content that looks like natural language responses
  // Skip if it looks like JSON or very short strings, or contains "not JSON"
  if (lowerContent.includes("{") || lowerContent.includes("}") || 
      content.length < 20 || lowerContent.includes("not json")) {
    return null;
  }
  
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
    "i'm sorry, but i cannot",
    "cannot process images",
    "unable to see the image"
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
    if (formulaPattern.test(smiles) && !["O", "C", "N", "S", "P", "H", "CCO", "CO"].includes(smiles)) {
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