// test/integration/molecular-accuracy.test.js - Real-world molecular accuracy testing
// Tests actual AtomPredictor responses against known chemical compositions
// Run with: npm test -- test/integration/molecular-accuracy.test.js

const Structuralizer = require("../../../backend/services/Structuralizer");

// Known molecular compositions for testing prompt accuracy
const REFERENCE_MATERIALS = {
  // Basic chemicals - core accuracy test
  basics: {
    "water": {
      expectedSMILES: ["O"],
      requiredNames: ["water"],
      forbiddenSMILES: ["H2O", "HOH"],
      accuracy_threshold: 95
    },
    "ethanol": {
      expectedSMILES: ["CCO"],
      requiredNames: ["ethanol"],
      forbiddenSMILES: ["C2H6O", "C2H5OH"],
      accuracy_threshold: 95
    },
    "sodium chloride": {
      expectedSMILES: ["[Na+].[Cl-]", "[Na+]", "[Cl-]"],
      requiredNames: ["sodium", "chloride"],
      forbiddenSMILES: ["NaCl"],
      accuracy_threshold: 85
    }
  },

  // Common beverages - realistic composition testing
  beverages: {
    "red wine": {
      expectedSMILES: ["CCO", "O", "OC(C(O)C(O)=O)C(O)=O"],
      requiredNames: ["ethanol", "water", "tartaric"],
      mayContain: ["glucose", "fructose", "malic", "tannin"],
      minComponents: 3,
      accuracy_threshold: 75
    },
    "black coffee": {
      expectedSMILES: ["O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"],
      requiredNames: ["water", "caffeine"],
      mayContain: ["chlorogenic", "acid"],
      minComponents: 2,
      accuracy_threshold: 80
    }
  },

  // Materials - complex composition handling
  materials: {
    "PET plastic bottle": {
      expectedSMILES: ["O=C(C1=CC=C(CO)C=C1)OC"],
      requiredNames: ["polyethylene", "terephthalate", "PET"],
      maxSMILESLength: 150,
      accuracy_threshold: 70
    },
    "limestone rock": {
      expectedSMILES: ["[Ca+2].[O-]C([O-])=O", "[Ca+2]", "[CO3-2]"],
      requiredNames: ["calcium", "carbonate"],
      forbiddenSMILES: ["CaCO3"],
      accuracy_threshold: 75
    }
  },

  // Biological materials - challenging but realistic
  biological: {
    "fresh apple": {
      expectedSMILES: ["O", "C(C(C(C(C(C=O)O)O)O)O)O"],
      requiredNames: ["water", "fructose", "glucose"],
      mayContain: ["malic", "cellulose", "pectin"],
      minComponents: 4,
      accuracy_threshold: 65
    },
    "olives": {
      // Expected SMILES moved from prompt example into test expectations
      expectedSMILES: [
        "OC1=CC=C(C=C1CCO)O", // Hydroxytyrosol
        "C1=CC(=CC=C1CCO)O",   // Tyrosol
        "CC(=CC1=CC(=CC=C1)O)CC=O", // Oleocanthal
        "C1=CC(=O)C2=C(O1)C(=CC(=C2O)O)C3=CC(=C(C=C3)O)O", // Luteolin
        "C1=CC(=O)C2=C(O1)C(=CC(=C2O)O)C3=CC=CC=C3", // Apigenin
        "CC1=CC(=C(C=C1OC2C(C(C(C(O2)CO)O)O)OC(=O)C=CC3=CC(=C(C=C3)O)O)O)OCC4C(C(C(C(O4)O)O)O)O" // Verbascoside (Acteoside)
      ],
      mayContain: ["oleuropein"],
      minComponents: 3,
      accuracy_threshold: 60
    }
  }
};

// Molecular analysis validator
class MolecularValidator {
  
  static validateSMILES(smiles) {
    if (!smiles || typeof smiles !== 'string') return false;
    
    // Must not be chemical formulas
    const formulaPattern = /^[A-Z][a-z]?\d*([A-Z][a-z]?\d*)*$/;
    if (formulaPattern.test(smiles) && !["O", "C", "N", "S", "P", "H"].includes(smiles)) {
      return { valid: false, error: "Chemical formula instead of SMILES" };
    }
    
    // Basic SMILES character validation
    const smilesPattern = /^[A-Za-z0-9\[\]()=#+\-\.@:\/\\%]+$/;
    if (!smilesPattern.test(smiles)) {
      return { valid: false, error: "Invalid SMILES characters" };
    }
    
    // Length validation
    if (smiles.length > 300) {
      return { valid: false, error: "SMILES too long" };
    }
    
    return { valid: true };
  }

  static analyzeMolecularAccuracy(result, reference) {
    const analysis = {
      score: 0,
      maxScore: 100,
      breakdown: {},
      issues: [],
      chemicals: result.chemicals || []
    };

    // Test 1: Expected SMILES present (40 points)
    if (reference.expectedSMILES) {
      const found = reference.expectedSMILES.some(expected => 
        result.chemicals.some(chemical => 
          chemical.smiles === expected || chemical.smiles.includes(expected)
        )
      );
      if (found) {
        analysis.score += 40;
        analysis.breakdown.expectedSMILES = "✓ Found";
      } else {
        analysis.issues.push(`Missing expected SMILES: ${reference.expectedSMILES.join(', ')}`);
        analysis.breakdown.expectedSMILES = "✗ Missing";
      }
    }

    // Test 2: Required names present (25 points) 
    if (reference.requiredNames) {
      const foundNames = reference.requiredNames.filter(name =>
        result.chemicals.some(chemical =>
          chemical.name.toLowerCase().includes(name.toLowerCase())
        )
      );
      const nameScore = (foundNames.length / reference.requiredNames.length) * 25;
      analysis.score += nameScore;
      analysis.breakdown.requiredNames = `${foundNames.length}/${reference.requiredNames.length} found`;
      
      if (foundNames.length < reference.requiredNames.length) {
        const missing = reference.requiredNames.filter(name => !foundNames.includes(name));
        analysis.issues.push(`Missing required names: ${missing.join(', ')}`);
      }
    }

    // Test 3: No forbidden SMILES (20 points)
    if (reference.forbiddenSMILES) {
      const hasForbidden = reference.forbiddenSMILES.some(forbidden =>
        result.chemicals.some(chemical => chemical.smiles === forbidden)
      );
      if (!hasForbidden) {
        analysis.score += 20;
        analysis.breakdown.forbiddenSMILES = "✓ None found";
      } else {
        analysis.issues.push(`Contains forbidden SMILES: ${reference.forbiddenSMILES.join(', ')}`);
        analysis.breakdown.forbiddenSMILES = "✗ Found forbidden";
      }
    } else {
      analysis.score += 20; // No forbidden SMILES to check
      analysis.breakdown.forbiddenSMILES = "N/A";
    }

    // Test 4: Minimum components (10 points)
    if (reference.minComponents) {
      if (result.chemicals.length >= reference.minComponents) {
        analysis.score += 10;
        analysis.breakdown.componentCount = `✓ ${result.chemicals.length} >= ${reference.minComponents}`;
      } else {
        analysis.issues.push(`Too few components: ${result.chemicals.length} < ${reference.minComponents}`);
        analysis.breakdown.componentCount = `✗ ${result.chemicals.length} < ${reference.minComponents}`;
      }
    } else {
      analysis.score += 10; // No minimum to check
      analysis.breakdown.componentCount = `${result.chemicals.length} components`;
    }

    // Test 5: SMILES validity (5 points)
    const invalidSMILES = result.chemicals.filter(chemical => {
      const validation = this.validateSMILES(chemical.smiles);
      return !validation.valid;
    });
    
    if (invalidSMILES.length === 0) {
      analysis.score += 5;
      analysis.breakdown.smilesValidity = "✓ All valid";
    } else {
      analysis.issues.push(`Invalid SMILES: ${invalidSMILES.map(c => c.smiles).join(', ')}`);
      analysis.breakdown.smilesValidity = `✗ ${invalidSMILES.length} invalid`;
    }

    analysis.score = Math.round(analysis.score);
    analysis.passed = analysis.score >= (reference.accuracy_threshold || 70);
    
    return analysis;
  }
}

// Skip these tests if no OpenAI API key is available
const skipIfNoApiKey = () => {
  if (!process.env.OPENAI_API_KEY) {
    console.log("⚠️  Skipping molecular accuracy tests - no OPENAI_API_KEY set");
    return true;
  }
  return false;
};

describe("Molecular Accuracy Integration Tests", () => {
  let atomPredictor;
  
  beforeAll(() => {
    if (skipIfNoApiKey()) return;
    atomPredictor = new Structuralizer(process.env.OPENAI_API_KEY);
  });

  describe("Basic Chemical Accuracy", () => {
    if (skipIfNoApiKey()) {
      test.skip("Skipping - no API key", () => {});
      return;
    }

    test.each(Object.entries(REFERENCE_MATERIALS.basics))(
      "should accurately analyze %s",
      async (material, reference) => {
        const result = await atomPredictor.structuralizeText(material);
        const analysis = MolecularValidator.analyzeMolecularAccuracy(result, reference);
        
        console.log(`\n=== ${material.toUpperCase()} ===`);
        console.log(`Score: ${analysis.score}% (threshold: ${reference.accuracy_threshold}%)`);
        console.log(`Breakdown:`, analysis.breakdown);
        if (analysis.issues.length > 0) {
          console.log(`Issues:`, analysis.issues);
        }
        console.log(`Chemicals found:`, result.chemicals.map(c => `${c.name}: ${c.smiles}`));

        expect(analysis.score).toBeGreaterThanOrEqual(reference.accuracy_threshold);
        expect(result.chemicals.length).toBeGreaterThan(0);
      },
      30000 // 30 second timeout for API calls
    );
  });

  describe("Beverage Composition Analysis", () => {
    if (skipIfNoApiKey()) {
      test.skip("Skipping - no API key", () => {});
      return;
    }

    test.each(Object.entries(REFERENCE_MATERIALS.beverages))(
      "should provide realistic composition for %s",
      async (beverage, reference) => {
        const result = await atomPredictor.structuralizeText(beverage);
        const analysis = MolecularValidator.analyzeMolecularAccuracy(result, reference);
        
        console.log(`\n=== ${beverage.toUpperCase()} ===`);
        console.log(`Score: ${analysis.score}% (threshold: ${reference.accuracy_threshold}%)`);
        console.log(`Components: ${result.chemicals.length}`);
        console.log(`Chemicals:`, result.chemicals.map(c => `${c.name}: ${c.smiles}`));
        
        expect(analysis.score).toBeGreaterThanOrEqual(reference.accuracy_threshold);
        expect(result.chemicals.length).toBeGreaterThanOrEqual(reference.minComponents || 1);
      },
      30000
    );
  });

  describe("Material Composition Analysis", () => {
    if (skipIfNoApiKey()) {
      test.skip("Skipping - no API key", () => {});
      return;
    }

    test.each(Object.entries(REFERENCE_MATERIALS.materials))(
      "should handle complex material %s",
      async (material, reference) => {
        const result = await atomPredictor.structuralizeText(material);
        const analysis = MolecularValidator.analyzeMolecularAccuracy(result, reference);
        
        console.log(`\n=== ${material.toUpperCase()} ===`);
        console.log(`Score: ${analysis.score}% (threshold: ${reference.accuracy_threshold}%)`);
        console.log(`Chemicals:`, result.chemicals.map(c => `${c.name}: ${c.smiles}`));
        
        // Check SMILES length constraints if specified
        if (reference.maxSMILESLength) {
          result.chemicals.forEach(chemical => {
            expect(chemical.smiles.length).toBeLessThanOrEqual(reference.maxSMILESLength);
          });
        }
        
        expect(analysis.score).toBeGreaterThanOrEqual(reference.accuracy_threshold);
      },
      30000
    );
  });

  describe("Biological Material Analysis", () => {
    if (skipIfNoApiKey()) {
      test.skip("Skipping - no API key", () => {});
      return;
    }

    test.each(Object.entries(REFERENCE_MATERIALS.biological))(
      "should provide realistic biological composition for %s",
      async (biological, reference) => {
        const result = await atomPredictor.structuralizeText(biological);
        const analysis = MolecularValidator.analyzeMolecularAccuracy(result, reference);
        
        console.log(`\n=== ${biological.toUpperCase()} ===`);
        console.log(`Score: ${analysis.score}% (threshold: ${reference.accuracy_threshold}%)`);
        console.log(`Components: ${result.chemicals.length}`);
        console.log(`Chemicals:`, result.chemicals.map(c => `${c.name}: ${c.smiles}`));
        
        // Biological materials should always contain water
        const hasWater = result.chemicals.some(c => 
          c.smiles === "O" || c.name.toLowerCase().includes("water")
        );
        expect(hasWater).toBe(true);
        
        expect(analysis.score).toBeGreaterThanOrEqual(reference.accuracy_threshold);
        expect(result.chemicals.length).toBeGreaterThanOrEqual(reference.minComponents || 1);
      },
      30000
    );
  });

  describe("SMILES Quality Validation", () => {
    if (skipIfNoApiKey()) {
      test.skip("Skipping - no API key", () => {});
      return;
    }

    test("should generate valid SMILES for all test materials", async () => {
      const testMaterials = ["water", "ethanol", "glucose", "caffeine"];
      let totalInvalid = 0;
      
      for (const material of testMaterials) {
        const result = await atomPredictor.analyzeText(material);
        
        result.chemicals.forEach(chemical => {
          const validation = MolecularValidator.validateSMILES(chemical.smiles);
          if (!validation.valid) {
            console.error(`Invalid SMILES for ${material}: ${chemical.smiles} - ${validation.error}`);
            totalInvalid++;
          }
        });
      }
      
      expect(totalInvalid).toBe(0);
    }, 60000);

    test("should not generate chemical formulas instead of SMILES", async () => {
      const result = await atomPredictor.structuralizeText("water and ethanol mixture");
      
      const forbiddenFormulas = ["H2O", "C2H6O", "C2H5OH"];
      const foundForbidden = result.chemicals.filter(chemical =>
        forbiddenFormulas.includes(chemical.smiles)
      );
      
      expect(foundForbidden).toHaveLength(0);
    }, 30000);
  });

  describe("Prompt Engineering Effectiveness", () => {
    if (skipIfNoApiKey()) {
      test.skip("Skipping - no API key", () => {});
      return;
    }

    test("should show improved accuracy with new prompt system", async () => {
      const testSuite = [
        ...Object.entries(REFERENCE_MATERIALS.basics),
        ...Object.entries(REFERENCE_MATERIALS.beverages)
      ];

      const results = [];
      
      for (const [material, reference] of testSuite) {
        const result = await atomPredictor.structuralizeText(material);
        const analysis = MolecularValidator.analyzeMolecularAccuracy(result, reference);
        
        results.push({
          material,
          score: analysis.score,
          passed: analysis.passed,
          threshold: reference.accuracy_threshold,
          chemicalCount: result.chemicals.length
        });
      }

      // Calculate metrics
      const averageScore = results.reduce((sum, r) => sum + r.score, 0) / results.length;
      const passRate = results.filter(r => r.passed).length / results.length;
      const averageChemicals = results.reduce((sum, r) => sum + r.chemicalCount, 0) / results.length;
      
      console.log("\n=== PROMPT ENGINEERING EFFECTIVENESS ===");
      console.log(`Average accuracy: ${averageScore.toFixed(1)}%`);
      console.log(`Pass rate: ${(passRate * 100).toFixed(1)}%`);
      console.log(`Average chemicals per analysis: ${averageChemicals.toFixed(1)}`);
      console.log(`Materials tested: ${results.length}`);
      
      // Individual results
      results.forEach(r => {
        const status = r.passed ? "✓" : "✗";
        console.log(`${status} ${r.material}: ${r.score}% (need ${r.threshold}%)`);
      });

      // Expectations based on our improvements
      expect(averageScore).toBeGreaterThanOrEqual(70); // Should be good
      expect(passRate).toBeGreaterThanOrEqual(0.6); // 60% should pass thresholds
      expect(averageChemicals).toBeGreaterThan(1); // Should provide multiple chemicals
      
    }, 120000); // 2 minute timeout for full suite
  });
});

// Export utilities for use in other tests
module.exports = {
  REFERENCE_MATERIALS,
  MolecularValidator
}; 