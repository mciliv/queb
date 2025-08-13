// Chemical Analysis Instructions for LLM
// This file contains proven prompting techniques extracted from git history

/**
 * Core Instruction Builder - Combines best practices from multiple iterations
 * 
 * KEY INSIGHTS FROM GIT HISTORY:
 * - Detailed examples significantly improve SMILES accuracy (from commit f1ddb4f)
 * - Explicit constraints prevent overly complex molecules (from commit 0a8f7b6) 
 * - Realistic chemical breakdowns work better than generic ones (from commit 0fe79e8)
 * - Comprehensive coverage with amounts adds context (from commit 40a48f2)
 */

function buildChemicalAnalysisInstructions() {
  return `Analyze object chemical composition. Return JSON.

ONLY physical objects: food, drinks, materials, minerals, chemicals.
REJECT: emotions, actions, nonsense text.

SMILES rules:
- Valid notation only: "[Na+].[Cl-]" for ionic compounds
- Standard database forms (PubChem, ChEBI)
- Proper syntax: atoms, bonds, rings, charges

Response format:
{
  "object": "Object name",
  "chemicals": [
    {"name": "Chemical name", "smiles": "SMILES notation"},
    {"name": "Chemical name", "smiles": "SMILES notation"}
  ]
}

Examples:
Water: "O", Ethanol: "CCO", NaCl: "[Na+].[Cl-]"
Glucose: "C(C(C(C(C(C=O)O)O)O)O)O"
Caffeine: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

Wine: Ethanol "CCO", Water "O", Tartaric acid "OC(C(O)C(O)=O)C(O)=O"

Avoid: H2O format, invalid syntax, unrealistic molecules`;
}

// Export for use in AtomPredictor
module.exports = { buildChemicalAnalysisInstructions }; 