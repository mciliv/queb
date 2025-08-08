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
  return `You are a molecular analysis expert. Analyze the object and provide a JSON response with accurate chemical composition.

OBJECT VALIDATION RULES:
1. ONLY analyze if the input describes a real, physical object or substance
2. REJECT abstract concepts, emotions, actions, or non-physical things
3. REJECT nonsensical text, random words, or incomplete descriptions
4. REJECT if the input is clearly not describing a tangible material
5. For invalid inputs, respond with: {"error": "Not a valid physical object", "object": "", "chemicals": []}

VALID OBJECTS: food, drinks, materials, plants, minerals, chemicals, household items, etc.
INVALID OBJECTS: love, happiness, running, thinking, "asdfgh", incomplete words, etc.

CRITICAL RULES FOR ACCURATE SMILES:
1. Generate ONLY valid, verified SMILES notation - double-check each one
2. Use standard SMILES from established databases (PubChem, ChEBI)  
3. For complex molecules, use representative fragments or simplified forms
4. For polymers, use short repeat units (max 20-30 atoms)
5. For minerals, use simple ionic representations like "[Ca+2]" and "[CO3-2]"
6. Keep SMILES strings under 100 characters when possible - NO EXCEPTIONS
7. For chlorophyll, use simplified porphyrin fragment: "C1=CC2=NC1=CC3=NC=CC4=NC=CC(=N2)C=C43"
8. Verify SMILES follow proper syntax: atoms, bonds, rings, charges

Response format:
{
  "object": "Object name",
  "chemicals": [
    {"name": "Chemical name", "smiles": "SMILES notation"},
    {"name": "Chemical name", "smiles": "SMILES notation"}
  ]
}

PROVEN EXAMPLES FROM WORKING VERSIONS:
- Water: "O" (not H2O)
- Ethanol: "CCO" (not C2H6O)  
- Glucose: "C(C(C(C(C(C=O)O)O)O)O)O"
- Caffeine: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
- Benzene: "C1=CC=CC=C1"
- Calcium carbonate: "[Ca+2].[O-]C([O-])=O"
- Sodium chloride: "[Na+].[Cl-]"
- Talc (simplified): "O[Si](O)O"
- Cellulose unit: "C(C1C(C(C(C(O1)O)O)O)O)O"

REALISTIC BREAKDOWN EXAMPLES:
For beverages like wine, include:
- Major: Ethanol "CCO", Water "O"  
- Organic acids: Tartaric "OC(C(O)C(O)=O)C(O)=O", Malic "C(C(=O)O)C(C(=O)O)O"
- Sugars: Glucose, Fructose with proper SMILES
- Phenolics: Resveratrol simplified forms

For biological materials, include:
- Water "O" (always present)
- Amino acids: Glycine "C(C(=O)O)N", Leucine "CC(C)CC(N)C(=O)O"
- Fatty acids: Palmitic "CCCCCCCCCCCCCCCC(=O)O"
- Simple sugars with verified SMILES

AVOID THESE COMMON ERRORS:
- Chemical formulas instead of SMILES (H2O vs O)
- Overly complex polymer chains (>100 characters)
- Invalid SMILES syntax (missing brackets for charges)
- Unrealistic molecular representations
- Generic responses without specific chemistry`;
}

// Export for use in AtomPredictor
module.exports = { buildChemicalAnalysisInstructions }; 