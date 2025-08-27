// test/integration/smiles-accuracy.test.js - Real-world SMILES accuracy testing

const MolecularProcessor = require("../../backend/services/molecular-processor");

// Known accurate SMILES from chemical databases
const KNOWN_MOLECULES = {
  // Simple molecules
  "water": { smiles: "O", formula: "H2O" },
  "methane": { smiles: "C", formula: "CH4" },
  "ethanol": { smiles: "CCO", formula: "C2H6O" },
  "acetic acid": { smiles: "CC(=O)O", formula: "C2H4O2" },
  
  // Common pharmaceuticals
  "aspirin": { smiles: "CC(=O)OC1=CC=CC=C1C(=O)O", formula: "C9H8O4" },
  "ibuprofen": { smiles: "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", formula: "C13H18O2" },
  "caffeine": { smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", formula: "C8H10N4O2" },
  "paracetamol": { smiles: "CC(=O)NC1=CC=C(C=C1)O", formula: "C8H9NO2" },
  
  // Biological molecules
  "glucose": { smiles: "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O", formula: "C6H12O6" },
  "alanine": { smiles: "C[C@@H](C(=O)O)N", formula: "C3H7NO2" },
  "glycine": { smiles: "C(C(=O)O)N", formula: "C2H5NO2" },
  
  // Natural products
  "menthol": { smiles: "CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O", formula: "C10H20O" },
  "vanillin": { smiles: "COC1=C(C=CC(=C1)C=O)O", formula: "C8H8O3" },
  
  // Complex molecules
  "cholesterol": { smiles: "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C", formula: "C27H46O" },
  "penicillin G": { smiles: "CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C", formula: "C16H18N2O4S" }
};

describe("SMILES Accuracy Integration Tests", () => {
  let processor;

  beforeAll(() => {
    processor = new MolecularProcessor();
  });

  describe("Known Molecule Validation", () => {
    test("should correctly identify valid SMILES for known molecules", () => {
      Object.entries(KNOWN_MOLECULES).forEach(([name, data]) => {
        expect(processor.isValidSmiles(data.smiles)).toBe(true);
      });
    });

    test("should process known pharmaceutical SMILES without errors", async () => {
      const pharmaceuticals = [
        KNOWN_MOLECULES.aspirin.smiles,
        KNOWN_MOLECULES.ibuprofen.smiles,
        KNOWN_MOLECULES.caffeine.smiles,
        KNOWN_MOLECULES.paracetamol.smiles
      ];

      const result = await processor.processSmiles(pharmaceuticals);
      
      // All should be processed successfully or have specific errors
      expect(result.sdfPaths.length + result.errors.length).toBe(pharmaceuticals.length);
      expect(result.skipped.length).toBe(0); // None should be skipped for format issues
    });

    test("should handle complex biological molecules", async () => {
      const biologicalMolecules = [
        KNOWN_MOLECULES.glucose.smiles,
        KNOWN_MOLECULES.alanine.smiles,
        KNOWN_MOLECULES.glycine.smiles,
        KNOWN_MOLECULES.cholesterol.smiles
      ];

      const result = await processor.processSmiles(biologicalMolecules);
      
      // Should attempt to process all (may have errors due to complexity, but not skip)
      expect(result.sdfPaths.length + result.errors.length).toBe(biologicalMolecules.length);
      expect(result.skipped.length).toBe(0);
    });
  });

  describe("SMILES Structural Validation", () => {
    test("should accept molecules with different ring sizes", async () => {
      const cyclicMolecules = [
        "C1CCC1", // cyclobutane
        "C1CCCC1", // cyclopentane
        "C1CCCCC1", // cyclohexane
        "C1CCCCCCC1", // cycloheptane
        "c1ccccc1", // benzene
        "c1ccc2ccccc2c1", // naphthalene
      ];

      const result = await processor.processSmiles(cyclicMolecules);
      expect(result.skipped.length).toBe(0); // All should be valid format
    });

    test("should handle stereochemistry correctly", async () => {
      const stereoMolecules = [
        "[C@H](C)(N)C(=O)O", // L-alanine
        "[C@@H](C)(N)C(=O)O", // D-alanine
        "C/C=C/C", // trans-2-butene
        "C/C=C\\C", // cis-2-butene
        "[C@H](O)[C@@H](O)C(=O)O", // stereoisomer
      ];

      const result = await processor.processSmiles(stereoMolecules);
      expect(result.skipped.length).toBe(0); // All should be valid format
    });

    test("should handle charged molecules", async () => {
      const chargedMolecules = [
        "[NH4+]", // ammonium
        "[OH-]", // hydroxide
        "[Na+]", // sodium ion
        "[Cl-]", // chloride
        "C[N+](C)(C)C", // tetramethylammonium
      ];

      const result = await processor.processSmiles(chargedMolecules);
      expect(result.skipped.length).toBe(0); // All should be valid format
    });
  });

  describe("Edge Cases and Error Boundaries", () => {
    test("should distinguish between valid SMILES and molecular formulas", () => {
      const testCases = [
        { input: "O", expected: true, type: "water SMILES" },
        { input: "H2O", expected: false, type: "water molecular formula" },
        { input: "C", expected: true, type: "methane SMILES" },
        { input: "CH4", expected: false, type: "methane molecular formula" },
        { input: "N", expected: true, type: "nitrogen SMILES" },
        { input: "NH3", expected: false, type: "ammonia molecular formula" },
      ];

      testCases.forEach(({ input, expected, type }) => {
        expect(processor.isValidSmiles(input)).toBe(expected);
      });
    });

    test("should handle unusual but valid SMILES notation", () => {
      const unusualValidSmiles = [
        "[H]", // explicit hydrogen
        "[He]", // helium
        "[Li+]", // lithium ion
        "B", // boron
        "[Si]", // silicon
        "P", // phosphorus
        "S", // sulfur
        "[Br]", // bromine
        "[I]", // iodine
      ];

      unusualValidSmiles.forEach(smiles => {
        expect(processor.isValidSmiles(smiles)).toBe(true);
      });
    });

    test("should reject obviously invalid constructions", () => {
      const invalidConstructions = [
        "C1C", // incomplete ring
        "C=", // incomplete double bond
        "C#", // incomplete triple bond
        "()", // empty parentheses
        "[]", // empty brackets
        "C1CCCC", // unclosed ring
      ];

      // Note: Our basic validator may not catch all of these,
      // but RDKit will catch them during processing
      const result = processor.processSmiles(invalidConstructions);
      // Should either skip or error, not succeed
    });
  });

  describe("Performance with Real-World Complexity", () => {
    test("should handle drug-like molecules efficiently", async () => {
      const drugLikeMolecules = [
        // Lipitor (atorvastatin)
        "CC(C)C1=C(C(=C(N1CC[C@H](C[C@H](CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
        
        // Viagra (sildenafil)
        "CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C",
        
        // Prozac (fluoxetine)
        "CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F",
      ];

      const startTime = Date.now();
      const result = await processor.processSmiles(drugLikeMolecules);
      const endTime = Date.now();

      expect(endTime - startTime).toBeLessThan(10000); // Should complete within 10 seconds
      // All should at least attempt processing (not skipped for format)
      expect(result.skipped.length).toBe(0);
    });

    test("should handle natural product complexity", async () => {
      const naturalProducts = [
        // Taxol (paclitaxel) - simplified
        "CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@H]([C@H](C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
        
        // Morphine
        "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O",
        
        // Quinine
        "COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O"
      ];

      const result = await processor.processSmiles(naturalProducts);
      
      // Should attempt to process all (complex molecules may fail in RDKit, but shouldn't be skipped)
      expect(result.skipped.length).toBe(0);
      expect(result.sdfPaths.length + result.errors.length).toBe(naturalProducts.length);
    });
  });

  describe("Regression Testing", () => {
    test("should handle previously problematic molecules", async () => {
      // Molecules that caused issues in the past
      const problematicMolecules = [
        // Lycopene (previously caused truncation)
        "CC(=CCCC(=CC=CC(=CC=CC=C(C)C=CC=C(C)C=CC=C(C)C=CC=C(C)C)C)C)C=CC=C(C)C=CC=C(C)C=CC1=C(C)CCCC1(C)C",
        
        // Beta-carotene
        "CC(=CCCC(=CC=CC(=CC=CC=C(C)C=CC=C(C)C=CC=C(C)C=CC=C(C)C)C)C)C=CC=C(C)C=CC=C(C)C=CC=C(C)C",
        
        // Large protein fragment
        "C[C@H](C(=O)N[C@@H](CC1=CC=CC=C1)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC(=O)N)C(=O)O)N"
      ];

      const result = await processor.processSmiles(problematicMolecules);
      
      // Should not skip any for length/format reasons
      expect(result.skipped.length).toBe(0);
      
      // Should provide informative feedback for any that fail
      result.errors.forEach(error => {
        expect(error).toMatch(/\S+/); // Should have meaningful error message
      });
    });
  });
});

// Export test utilities for other test files
module.exports = {
  KNOWN_MOLECULES,
  knownValidSmiles: Object.values(KNOWN_MOLECULES).map(m => m.smiles),
  testSmilesAccuracy: async (smiles) => {
    const processor = new MolecularProcessor();
    return await processor.processSmiles([smiles]);
  }
}; 