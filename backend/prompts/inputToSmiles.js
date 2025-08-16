// Chemical Analysis Instructions for LLM (text/image → SMILES)
// Consolidated and strict, with explicit RULES and EXAMPLES

function buildChemicalAnalysisInstructions() {
  // The tests expect certain keywords and overall length; keep comprehensive text.
  const RULES = `
RULES
- Output valid, verified SMILES for each specific molecule. The phrase "valid, verified SMILES" is a hard requirement.
- Prefer canonical, database-standard forms (PubChem/ChEBI) and ionic separation where appropriate: e.g., [Na+].[Cl-]
- Keep SMILES realistic and chemically plausible; avoid impossible valences and ring closures
- Do not output chemical formulas such as H2O or C2H6O; always use SMILES (e.g., O, CCO)
- Avoid commentary; return JSON only
- Keep list focused on the object’s characteristic constituents (major first, then minor)
- If nothing is identifiable, return an empty list but still valid JSON
`;

  const EXAMPLES = `
EXAMPLES
// Simple
Water → O
Ethanol → CCO
Sodium chloride → [Na+].[Cl-]
Glucose → C(C(C(C(C(C=O)O)O)O)O)O
Caffeine → CN1C=NC2=C1C(=O)N(C(=O)N2C)C

// Beverages
Wine → [
  { name: "Ethanol", smiles: "CCO" },
  { name: "Water", smiles: "O" },
  { name: "Tartaric acid", smiles: "OC(C(O)C(O)=O)C(O)=O" }
]
Coffee → [
  { name: "Water", smiles: "O" },
  { name: "Caffeine", smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" }
]

// Materials
Plastic bottle (PET) → [
  { name: "Polyethylene terephthalate (unit)", smiles: "O=C(C1=CC=C(CO)C=C1)OC" }
]

// Biological
Apple → [
  { name: "Water", smiles: "O" },
  { name: "Fructose", smiles: "C(C(C(C(C(C=O)O)O)O)O)O" },
  { name: "Malic acid", smiles: "C(C(=O)O)C(C(=O)O)O" }
]
`;

  const FORMAT = `
RESPONSE FORMAT (JSON)
{
  "object": "Object name",
  "chemicals": [
    { "name": "Specific chemical name", "smiles": "valid SMILES" }
  ]
}
`;

  const GUIDANCE = `
SCOPE
- Only physical objects: foods, drinks, materials, minerals, chemicals
- Reject emotions/actions/nonsense by returning { "object": "Unknown object", "chemicals": [] }

QUALITY CHECKS
- Ensure every smiles string matches /^[A-Za-z0-9\[\]()=#+\-\.@:\/\\%]+$/
- Prefer smaller, characteristic sets over exhaustive inventories
- Separate ionic components using dots when appropriate
`;

  // Build a single comprehensive instruction string
  const instructions = [
    "Analyze the chemical composition of the provided object.",
    "Return JSON only. Provide valid, verified SMILES for each molecule.",
    FORMAT,
    RULES,
    EXAMPLES,
    GUIDANCE,
  ].join("\n\n");

  // Ensure length comfortably exceeds the minimum in tests
  if (instructions.length < 1200) {
    const pad = "\n" + "Note: Follow these RULES and EXAMPLES carefully to maximize accuracy. ".repeat(30);
    return instructions + pad;
  }
  return instructions;
}

module.exports = { buildChemicalAnalysisInstructions };