// Minimal validation patterns - only block truly problematic inputs
export const VALIDATION_PATTERNS = {
  // Only block inputs that are clearly not useful
  onlySymbols: /^[^a-zA-Z0-9\s]+$/,  // Only symbols/punctuation
  emptyOrWhitespace: /^\s*$/          // Empty or only whitespace
};

// Positive patterns for likely valid inputs
export const POSITIVE_PATTERNS = {
  chemicalNames: /^(acid|alcohol|amine|benzene|methyl|ethyl|propyl|butyl|phenyl|hydroxyl|carboxyl).*$/i,
  biologicalTerms: /^(protein|enzyme|hormone|vitamin|amino|glucose|fructose|cellulose|lignin|starch).*$/i,
  foodItems: /^(apple|banana|orange|lemon|garlic|onion|tomato|potato|carrot|lettuce|spinach|broccoli|cheese|milk|butter|oil|sugar|salt|pepper).*$/i,
  materials: /^(plastic|metal|wood|glass|ceramic|rubber|silicon|carbon|diamond|graphite|steel|aluminum|copper|iron|gold|silver).*$/i,
  minerals: /^(quartz|feldspar|mica|calcite|gypsum|halite|pyrite|hematite|magnetite|olivine).*$/i
};

// Smart validation function
export const validateInput = (text) => {
  if (!text || !text.trim()) {
    return 'Enter a thing to structuralize';
  }
  
  const trimmed = text.trim().toLowerCase();
  
  if (trimmed.length < 2) {
    return 'Input must be at least 2 characters long';
  }
  if (trimmed.length > 500) {
    return 'Input must be less than 500 characters';
  }
  
  // Check negative patterns
  const invalidPatterns = Object.values(VALIDATION_PATTERNS);
  if (invalidPatterns.some(pattern => pattern.test(trimmed))) {
    return 'Please describe a real, physical object (food, materials, plants, etc.)';
  }
  
  // Additional smart checks
  const words = trimmed.split(/\s+/);
  if (words.length === 1 && words[0].length < 3 && !/^(oil|tea|wax|ice|air)$/.test(words[0])) {
    return 'Please provide a more specific description';
  }
  
  // Check for repeated characters (like "aaaa" or "123123")
  if (/^(..)\1{2,}$/.test(trimmed) || /^(.)\1{4,}$/.test(trimmed)) {
    return 'Please enter a valid substance name';
  }
  
  return null;
};

// Input suggestion scorer
export const scoreInputSuggestion = (suggestion, input) => {
  const lowerSuggestion = suggestion.toLowerCase();
  const lowerInput = input.toLowerCase();
  
  let score = 0;
  
  // Exact match gets highest score
  if (lowerSuggestion === lowerInput) {
    return 100;
  }
  
  // Starts with input gets high score
  if (lowerSuggestion.startsWith(lowerInput)) {
    score += 80;
  }
  
  // Contains input gets medium score
  if (lowerSuggestion.includes(lowerInput)) {
    score += 40;
  }
  
  // Positive pattern match gets bonus
  const positivePatterns = Object.values(POSITIVE_PATTERNS);
  if (positivePatterns.some(pattern => pattern.test(lowerSuggestion))) {
    score += 10;
  }
  
  // Length similarity bonus
  const lengthRatio = Math.min(lowerInput.length, lowerSuggestion.length) / 
                     Math.max(lowerInput.length, lowerSuggestion.length);
  score += lengthRatio * 10;
  
  return score;
};