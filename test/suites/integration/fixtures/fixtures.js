// test/fixtures/fixtures.js - Test fixtures and mock data
const fs = require("fs");
const path = require("path");

// ==================== TEST FIXTURES ====================

// Common test molecules with detailed information
const TEST_MOLECULES = {
  // Simple molecules
  water: {
    smiles: "O",
    name: "Water",
    formula: "H2O",
    description: "Universal solvent",
    molarMass: 18.015,
  },

  ethanol: {
    smiles: "CCO",
    name: "Ethanol",
    formula: "C2H6O",
    description: "Alcohol found in beverages",
    molarMass: 46.07,
  },

  methane: {
    smiles: "C",
    name: "Methane",
    formula: "CH4",
    description: "Simplest hydrocarbon",
    molarMass: 16.04,
  },

  // Complex molecules
  caffeine: {
    smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    name: "Caffeine",
    formula: "C8H10N4O2",
    description: "Stimulant found in coffee",
    molarMass: 194.19,
  },

  aspirin: {
    smiles: "CC(=O)OC1=CC=CC=C1C(=O)O",
    name: "Aspirin",
    formula: "C9H8O4",
    description: "Pain reliever and anti-inflammatory",
    molarMass: 180.16,
  },

  glucose: {
    smiles: "C(C1C(C(C(C(O1)O)O)O)O)O",
    name: "Glucose",
    formula: "C6H12O6",
    description: "Simple sugar",
    molarMass: 180.16,
  },

  // Organic compounds
  benzene: {
    smiles: "C1=CC=CC=C1",
    name: "Benzene",
    formula: "C6H6",
    description: "Aromatic hydrocarbon",
    molarMass: 78.11,
  },

  toluene: {
    smiles: "CC1=CC=CC=C1",
    name: "Toluene",
    formula: "C7H8",
    description: "Aromatic hydrocarbon solvent",
    molarMass: 92.14,
  },
};

// Test objects and their expected molecular components
const TEST_OBJECTS = {
  coffee: {
    object: "coffee",
    expectedMolecules: ["caffeine", "water", "chlorogenic acid"],
    smiles: [
      "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", // caffeine
      "O", // water
      "C1=CC(=C(C=C1C=CC(=O)O)O)O", // chlorogenic acid
    ],
    description: "Brewed coffee beverage",
  },

  wine: {
    object: "wine",
    expectedMolecules: ["ethanol", "water", "tartaric acid"],
    smiles: [
      "CCO", // ethanol
      "O", // water
      "OC(C(O)C(O)=O)C(O)=O", // tartaric acid
    ],
    description: "Fermented grape beverage",
  },

  aspirinTablet: {
    object: "aspirin tablet",
    expectedMolecules: ["aspirin", "cellulose", "starch"],
    smiles: [
      "CC(=O)OC1=CC=CC=C1C(=O)O", // aspirin
      "C1C(C(C(C(O1)CO)O)O)O", // cellulose
      "C(C1C(C(C(C(O1)O)O)O)O)O", // starch
    ],
    description: "Pharmaceutical tablet",
  },

  plasticBottle: {
    object: "plastic bottle",
    expectedMolecules: ["polyethylene terephthalate", "additives"],
    smiles: [
      "C(C(=O)OC1=CC=C(C=C1)C2=CC=C(C=C2)OC(=O)C)OC(=O)C3=CC=C(C=C3)C4=CC=C(C=C4)OC(=O)C", // PET
      "CCCCCCCCCCCCCCC(=O)O", // palmitic acid (additive)
    ],
    description: "Plastic container",
  },

  wood: {
    object: "wood",
    expectedMolecules: ["cellulose", "lignin", "hemicellulose"],
    smiles: [
      "C1C(C(C(C(O1)CO)O)O)O", // cellulose
      "COC1=CC(=CC(=C1O)OC)C=CC(=O)O", // lignin component
      "C(C1C(C(C(C(O1)O)O)O)O)O", // hemicellulose
    ],
    description: "Natural wood material",
  },
};

// Mock image data for testing
const MOCK_IMAGES = {
  // 1x1 pixel images for testing
  blackSquare: {
    base64:
      "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==",
    description: "1x1 black pixel",
  },

  whiteSquare: {
    base64:
      "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNk+M9QDwADhgGAWjR9awAAAABJRU5ErkJggg==",
    description: "1x1 white pixel",
  },

  // Larger test images
  colorGradient: {
    base64:
      "iVBORw0KGgoAAAANSUhEUgAAAAoAAAAKCAYAAACNMs+9AAAABUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==",
    description: "10x10 gradient for testing",
  },
};

// Test response templates
const RESPONSE_TEMPLATES = {
  successful: {
    imageMolecules: (object, smiles) => ({
      output: {
        object,
        smiles,
        _test: {
          fixture: true,
          timestamp: new Date().toISOString(),
        },
      },
    }),

    objectMolecules: (smiles) => ({
      output: {
        smiles,
        _test: {
          fixture: true,
          timestamp: new Date().toISOString(),
        },
      },
    }),

    generateSdfs: (sdfPaths) => ({
      message: "Files generated",
      sdfPaths,
      _test: {
        fixture: true,
        timestamp: new Date().toISOString(),
      },
    }),
  },

  error: {
    missingData: (field) => ({
      error: `${field} is required`,
      _test: {
        fixture: true,
        errorType: "validation",
        timestamp: new Date().toISOString(),
      },
    }),

    serverError: (message) => ({
      error: message,
      _test: {
        fixture: true,
        errorType: "server",
        timestamp: new Date().toISOString(),
      },
    }),
  },
};

// ==================== HELPER FUNCTIONS ====================

/**
 * Get a test molecule by name
 */
function getTestMolecule(name) {
  return TEST_MOLECULES[name];
}

/**
 * Get a test object by name
 */
function getTestObject(name) {
  return TEST_OBJECTS[name];
}

/**
 * Get random test SMILES
 */
function getRandomSmiles() {
  const molecules = Object.values(TEST_MOLECULES);
  return molecules[Math.floor(Math.random() * molecules.length)].smiles;
}

/**
 * Get random test object
 */
function getRandomObject() {
  const objects = Object.keys(TEST_OBJECTS);
  return objects[Math.floor(Math.random() * objects.length)];
}

/**
 * Create mock SDF file content
 */
function createMockSdf(smiles, moleculeName = "TestMolecule") {
  return `
  Mrv2114 ${new Date().toISOString().slice(0, 10).replace(/-/g, "")}

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <SMILES>
${smiles}

> <Name>
${moleculeName}

 > <_test_file>
 true

$$$$
`;
}

/**
 * Create test request data
 */
function createTestRequest(type, options = {}) {
  const templates = {
    imageMolecules: {
      imageBase64: options.imageBase64 || MOCK_IMAGES.blackSquare.base64,
      croppedImageBase64:
        options.croppedImageBase64 || MOCK_IMAGES.whiteSquare.base64,
      x: options.x || 100,
      y: options.y || 100,
      ...options,
    },

    objectMolecules: {
      object: options.object || getRandomObject(),
      ...options,
    },

    generateSdfs: {
      smiles: options.smiles || [getRandomSmiles()],
      overwrite: options.overwrite || false,
      ...options,
    },
  };

  return templates[type] || {};
}

/**
 * Validate test response
 */
function validateTestResponse(response, expectedType) {
  const validators = {
    imageMolecules: (res) => {
      return (
        res.output && res.output.object && Array.isArray(res.output.smiles)
      );
    },

    objectMolecules: (res) => {
      return res.output && Array.isArray(res.output.smiles);
    },

    generateSdfs: (res) => {
      return res.message && Array.isArray(res.sdfPaths);
    },
  };

  const validator = validators[expectedType];
  return validator ? validator(response) : false;
}

// ==================== EXPORTS ====================
module.exports = {
  TEST_MOLECULES,
  TEST_OBJECTS,
  MOCK_IMAGES,
  RESPONSE_TEMPLATES,
  getTestMolecule,
  getTestObject,
  getRandomSmiles,
  getRandomObject,
  createMockSdf,
  createTestRequest,
  validateTestResponse,
};
