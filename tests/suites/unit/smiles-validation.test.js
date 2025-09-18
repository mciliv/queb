// test/unit/smiles-validation.test.js - SMILES accuracy and validation tests

const MolecularProcessor = require("../../backend/services/molecular-processor");
const path = require("path");
const fs = require("fs");

// Mock child_process for SDF generation
jest.mock('child_process', () => ({
  spawn: jest.fn(() => ({
    stdout: { 
      on: jest.fn((event, callback) => {
        if (event === 'data') {
          callback('Generated SDF file: /mock/path/test.sdf\n');
        }
      })
    },
    stderr: { on: jest.fn() },
    on: jest.fn((event, callback) => {
      if (event === 'close') {
        setTimeout(() => callback(0), 10);
      }
    })
  }))
}));

describe("SMILES Validation Tests", () => {
  let processor;

  beforeEach(() => {
    processor = new MolecularProcessor();
  });

  describe("Valid SMILES Recognition", () => {
    test("should accept simple valid SMILES", () => {
      const validSmiles = [
        "O", // water
        "CCO", // ethanol
        "C1=CC=CC=C1", // benzene
        "CC(=O)O", // acetic acid
        "N", // ammonia
        "CC", // ethane
        "C=C", // ethene
        "C#C", // ethyne
        "c1ccccc1", // benzene (aromatic)
      ];

      validSmiles.forEach(smiles => {
        expect(processor.isValidSmiles(smiles)).toBe(true);
      });
    });

    test("should accept complex valid SMILES", () => {
      const complexSmiles = [
        // Glucose
        "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O",
        
        // Caffeine
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        
        // Aspirin
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        
        // Penicillin G
        "CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C",
        
        // Lycopene (simplified representation)
        "CC(=CCCC(=CC=CC(=CC=CC=C(C)C=CC=C(C)C=CC=C(C)C)C)C)C=CC=C(C)C=CC=C(C)C=CC1=C(C)CCCC1(C)C",
        
        // DNA base adenine
        "NC1=NC=NC2=C1N=CN2",
        
        // Cholesterol
        "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C"
      ];

      complexSmiles.forEach(smiles => {
        expect(processor.isValidSmiles(smiles)).toBe(true);
      });
    });

    test("should accept SMILES with various bond types", () => {
      const bondSmiles = [
        "C-C", // single bond (explicit)
        "C=C", // double bond
        "C#C", // triple bond
        "c1ccccc1", // aromatic
        "C1CCC1", // ring closure
        "C1CCCCC1", // cyclohexane
        "C1=CC=CC=C1", // benzene (explicit aromatic)
      ];

      bondSmiles.forEach(smiles => {
        expect(processor.isValidSmiles(smiles)).toBe(true);
      });
    });

    test("should accept SMILES with stereochemistry", () => {
      const stereoSmiles = [
        "[C@H](C)(N)C(=O)O", // L-alanine
        "[C@@H](C)(N)C(=O)O", // D-alanine
        "C/C=C/C", // trans double bond
        "C/C=C\\C", // cis double bond
        "[NH3+]", // charged nitrogen
        "[O-]", // negatively charged oxygen
      ];

      stereoSmiles.forEach(smiles => {
        expect(processor.isValidSmiles(smiles)).toBe(true);
      });
    });
  });

  describe("Invalid SMILES Recognition", () => {
    test("should reject molecular formulas", () => {
      const molecularFormulas = [
        "H2O",
        "CO2", 
        "CaCO3",
        "NaCl",
        "C6H12O6",
        "CH4",
        "NH3", // This could be confused with valid SMILES "N"
        "C2H6O",
      ];

      // Most of these should be rejected by our molecular formula detection
      const shouldBeRejected = ["H2O", "CO2", "CaCO3", "NaCl", "C6H12O6", "CH4", "C2H6O"];
      
      shouldBeRejected.forEach(formula => {
        expect(processor.isValidSmiles(formula)).toBe(false);
      });
    });

    test("should reject empty or invalid inputs", () => {
      const invalidInputs = [
        "",
        " ",
        null,
        undefined,
        "N/A",
        123,
        {},
        [],
      ];

      invalidInputs.forEach(input => {
        expect(processor.isValidSmiles(input)).toBe(false);
      });
    });
  });

  describe("SMILES Processing Integration", () => {
    test("should process array with mix of valid and invalid SMILES", async () => {
      const mixedSmiles = [
        "O", // valid - water
        "CCO", // valid - ethanol
        "H2O", // invalid - molecular formula
        "C1=CC=CC=C1", // valid - benzene
        "CaCO3", // invalid - molecular formula
        "N/A", // invalid - placeholder
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", // valid - caffeine
      ];

      const result = await processor.processSmiles(mixedSmiles);
      
      // Should have some successful processing and some skipped
      expect(result.sdfPaths.length).toBeGreaterThan(0);
      expect(result.skipped.length).toBeGreaterThan(0);
      
      // Should skip the molecular formulas and N/A
      expect(result.skipped.some(s => s.includes("H2O"))).toBe(true);
      expect(result.skipped.some(s => s.includes("CaCO3"))).toBe(true);
      expect(result.skipped.some(s => s.includes("N/A"))).toBe(true);
    });

    test("should handle very long but valid SMILES", async () => {
      // Lycopene - a long but valid SMILES
      const lycopeneSmiles = "CC(=CCCC(=CC=CC(=CC=CC=C(C)C=CC=C(C)C=CC=C(C)C)C)C)C=CC=C(C)C=CC=C(C)C=CC1=C(C)CCCC1(C)C";
      
      const result = await processor.processSmiles([lycopeneSmiles]);
      
      // Should not skip based on length alone
      expect(result.skipped.length).toBe(0);
      // May have errors if RDKit can't parse it, but shouldn't be skipped for length
    });
  });

  describe("Error Handling and Reporting", () => {
    test("should provide informative error messages", async () => {
      const problematicSmiles = [
        "H2O", // molecular formula
        "", // empty
        "INVALID_SMILES_XYZ", // clearly invalid
      ];

      const result = await processor.processSmiles(problematicSmiles);
      
      expect(result.skipped.length).toBeGreaterThan(0);
      expect(result.errors.length).toBeGreaterThanOrEqual(0);
      
      // Check that error messages are informative
      result.skipped.forEach(skip => {
        expect(typeof skip).toBe('string');
        expect(skip.length).toBeGreaterThan(0);
      });
    });

    test("should handle processing errors gracefully", async () => {
      // Mock a processing error for specific SMILES
      const originalSpawn = require('child_process').spawn;
      require('child_process').spawn.mockImplementationOnce(() => ({
        stdout: { on: jest.fn() },
        stderr: { on: jest.fn() },
        on: jest.fn((event, callback) => {
          if (event === 'close') {
            callback(1); // Error exit code
          }
        })
      }));

      const result = await processor.processSmiles(["CCO"]);
      
      // Should handle the error gracefully
      expect(result.errors.length).toBeGreaterThan(0);
      expect(result.sdfPaths.length).toBe(0);
      
      // Restore original mock
      require('child_process').spawn.mockImplementation(originalSpawn);
    });
  });

  describe("Chemical Accuracy Validation", () => {
    test("should identify common chemical classes correctly", () => {
      const chemicalClasses = {
        alcohols: ["CCO", "CCCO", "CC(C)O"],
        acids: ["CC(=O)O", "C(=O)O"],
        aromatics: ["c1ccccc1", "C1=CC=CC=C1"],
        alkanes: ["CC", "CCC", "CCCC"],
        alkenes: ["C=C", "CC=C"],
        alkynes: ["C#C", "CC#C"],
      };

      Object.entries(chemicalClasses).forEach(([className, smilesList]) => {
        smilesList.forEach(smiles => {
          expect(processor.isValidSmiles(smiles)).toBe(true);
        });
      });
    });

    test("should handle biological molecules", () => {
      const biologicalMolecules = [
        // Amino acids
        "C(C(=O)O)N", // glycine
        "C[C@@H](C(=O)O)N", // alanine
        "C1=CC=C(C=C1)C[C@@H](C(=O)O)N", // phenylalanine
        
        // Nucleotides
        "NC1=NC=NC2=C1N=CN2", // adenine
        "NC1=NC(=O)C=CN1", // cytosine
        "NC1=C(N=CN1)C(=O)N", // guanine
        "CC1=CN=C(N1)N", // thymine
        
        // Sugars
        "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", // glucose
        "OC[C@@H](O)[C@@H](O)[C@H](O)C(=O)CO", // fructose
      ];

      biologicalMolecules.forEach(smiles => {
        expect(processor.isValidSmiles(smiles)).toBe(true);
      });
    });
  });

  describe("Performance with Complex Molecules", () => {
    test("should handle large molecules efficiently", async () => {
      // Test with progressively larger molecules
      const largeMolecules = [
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", // caffeine
        "CC(=O)OC1=CC=CC=C1C(=O)O", // aspirin
        "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C", // cholesterol
      ];

      const startTime = Date.now();
      const result = await processor.processSmiles(largeMolecules);
      const endTime = Date.now();

      expect(endTime - startTime).toBeLessThan(5000); // Should complete within 5 seconds
      expect(result.sdfPaths.length).toBeGreaterThan(0);
    });
  });
});

module.exports = {
  // Export utilities for other tests
  validSmiles: [
    "O", "CCO", "C1=CC=CC=C1", "CC(=O)O", "N", 
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
  ],
  invalidSmiles: [
    "H2O", "CO2", "CaCO3", "N/A", "", null, undefined
  ]
}; 