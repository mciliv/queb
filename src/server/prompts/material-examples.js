// Proven Chemical Analysis Examples
// Curated from successful prompts in git history

/**
 * Comprehensive Examples Library
 * 
 * KEY INSIGHTS FROM GIT HISTORY:
 * - Detailed examples dramatically improve accuracy (commit f1ddb4f)
 * - Realistic breakdowns for complex materials work best (commit be3ecc0)  
 * - Specific examples guide LLM toward correct SMILES format (all commits)
 * - Showing both major and minor components improves completeness
 */

// Beverages - From highly successful prompts (commit f1ddb4f)
const BEVERAGE_EXAMPLES = {
  water: {
    object: "Tap water",
    chemicals: [
      {"name": "Water", "smiles": "O"},
      {"name": "Sodium ion", "smiles": "[Na+]"},
      {"name": "Chloride ion", "smiles": "[Cl-]"}, 
      {"name": "Calcium ion", "smiles": "[Ca+2]"},
      {"name": "Magnesium ion", "smiles": "[Mg+2]"}
    ]
  },
  
  wine: {
    object: "Wine",
    chemicals: [
      {"name": "Ethanol", "smiles": "CCO"},
      {"name": "Water", "smiles": "O"},
      {"name": "Tartaric acid", "smiles": "OC(C(O)C(O)=O)C(O)=O"},
      {"name": "Glucose", "smiles": "C(C(C(C(C(C=O)O)O)O)O)O"},
      {"name": "Fructose", "smiles": "C(C(C(C(C(C=O)O)O)O)O)O"},
      {"name": "Malic acid", "smiles": "C(C(=O)O)C(C(=O)O)O"}
    ]
  },

  coffee: {
    object: "Coffee",
    chemicals: [
      {"name": "Water", "smiles": "O"},
      {"name": "Caffeine", "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"},
      {"name": "Chlorogenic acid", "smiles": "C1=CC(=C(C=C1)O)O"},
      {"name": "Glucose", "smiles": "C(C(C(C(C(C=O)O)O)O)O)O"}
    ]
  }
};

// Materials and objects - Proven effective examples
const MATERIAL_EXAMPLES = {
  plastic: {
    object: "Plastic bottle (PET)",
    chemicals: [
      {"name": "Polyethylene terephthalate unit", "smiles": "O=C(C1=CC=C(CO)C=C1)OC"},
      {"name": "Ethylene glycol", "smiles": "OCCO"},
      {"name": "Terephthalic acid", "smiles": "O=C(O)C1=CC=C(C(=O)O)C=C1"}
    ]
  },

  limestone: {
    object: "Limestone (calcium carbonate)", 
    chemicals: [
      {"name": "Calcium carbonate", "smiles": "[Ca+2].[O-]C([O-])=O"},
      {"name": "Magnesium carbonate", "smiles": "[Mg+2].[O-]C([O-])=O"},
      {"name": "Silicon dioxide", "smiles": "O=[Si]=O"}
    ]
  },

  wood: {
    object: "Wood",
    chemicals: [
      {"name": "Cellulose unit", "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O"},
      {"name": "Lignin unit", "smiles": "C1=CC(=C(C=C1)O)O"}, 
      {"name": "Water", "smiles": "O"}
    ]
  }
};

// Biological materials - Realistic and proven effective
const BIOLOGICAL_EXAMPLES = {
  fruit: {
    object: "Apple",
    chemicals: [
      {"name": "Water", "smiles": "O"},
      {"name": "Fructose", "smiles": "C(C(C(C(C(C=O)O)O)O)O)O"},
      {"name": "Glucose", "smiles": "C(C(C(C(C(C=O)O)O)O)O)O"},
      {"name": "Malic acid", "smiles": "C(C(=O)O)C(C(=O)O)O"},
      {"name": "Cellulose unit", "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O"},
      {"name": "Pectin unit", "smiles": "CC(=O)OC1C(OC(CO)C(O)C1O)OC"}
    ]
  },

  meat: {
    object: "Cooked meat",
    chemicals: [
      {"name": "Water", "smiles": "O"},
      {"name": "Glycine", "smiles": "C(C(=O)O)N"},
      {"name": "Alanine", "smiles": "CC(C(=O)O)N"},
      {"name": "Leucine", "smiles": "CC(C)CC(C(=O)O)N"},
      {"name": "Palmitic acid", "smiles": "CCCCCCCCCCCCCCCC(=O)O"},
      {"name": "Oleic acid", "smiles": "CCCCCCCCC=CCCCCCCCC(=O)O"}
    ]
  }
};

// Simple chemicals - Critical for accuracy
const SIMPLE_CHEMICAL_EXAMPLES = {
  basic: [
    {"name": "Water", "smiles": "O"},
    {"name": "Ethanol", "smiles": "CCO"},
    {"name": "Methane", "smiles": "C"},
    {"name": "Benzene", "smiles": "C1=CC=CC=C1"},
    {"name": "Glucose", "smiles": "C(C(C(C(C(C=O)O)O)O)O)O"}
  ],
  
  salts: [
    {"name": "Sodium chloride", "smiles": "[Na+].[Cl-]"},
    {"name": "Calcium carbonate", "smiles": "[Ca+2].[O-]C([O-])=O"},
    {"name": "Potassium nitrate", "smiles": "[K+].[O-][N+]([O-])=O"},
    {"name": "Magnesium sulfate", "smiles": "[Mg+2].[O-]S([O-])(=O)=O"}
  ],

  acids: [
    {"name": "Acetic acid", "smiles": "CC(=O)O"},
    {"name": "Citric acid", "smiles": "C(C(=O)O)C(O)(C(=O)O)C(=O)O"},
    {"name": "Lactic acid", "smiles": "C(C(=O)O)O"},
    {"name": "Tartaric acid", "smiles": "OC(C(O)C(O)=O)C(O)=O"}
  ]
};

/**
 * Generate context-aware examples for the LLM prompt
 * This helps the LLM understand the expected format and accuracy level
 */
function getRelevantExamples(objectType = 'general') {
  const baseExamples = `
PROVEN ACCURATE EXAMPLES:

Basic molecules:
- Water: "O" (NOT H2O)
- Ethanol: "CCO" (NOT C2H6O)
- Glucose: "C(C(C(C(C(C=O)O)O)O)O)O"
- Caffeine: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
- Benzene: "C1=CC=CC=C1"

Ionic compounds:
- Sodium chloride: "[Na+].[Cl-]"
- Calcium carbonate: "[Ca+2].[O-]C([O-])=O"
- Magnesium sulfate: "[Mg+2].[O-]S([O-])(=O)=O"

Complex materials:
- Cellulose unit: "C(C1C(C(C(C(O1)O)O)O)O)O"
- Protein amino acid: "C(C(=O)O)N" (glycine)
- Fatty acid: "CCCCCCCCCCCCCCCC(=O)O" (palmitic)`;

  // Add specific examples based on context
  if (objectType.includes('beverage') || objectType.includes('drink')) {
    return baseExamples + `

Beverage-specific examples:
- Wine: Ethanol "CCO", Water "O", Tartaric acid "OC(C(O)C(O)=O)C(O)=O"
- Coffee: Water "O", Caffeine "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
- Beer: Ethanol "CCO", Water "O", Glucose sugars`;
  }

  if (objectType.includes('food') || objectType.includes('fruit')) {
    return baseExamples + `

Food-specific examples:
- Fruits: Water "O", Fructose sugars, Malic acid "C(C(=O)O)C(C(=O)O)O"
- Vegetables: Water "O", Cellulose "C(C1C(C(C(C(O1)O)O)O)O)O", various minerals`;
  }

  return baseExamples;
}

module.exports = {
  BEVERAGE_EXAMPLES,
  MATERIAL_EXAMPLES, 
  BIOLOGICAL_EXAMPLES,
  SIMPLE_CHEMICAL_EXAMPLES,
  getRelevantExamples
}; 