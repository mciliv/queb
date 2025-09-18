// Test utilities for integration tests
const fs = require('fs');
const path = require('path');

class TestFileManager {
  constructor() {
    this.tempFiles = [];
  }

  createTempFile(content, extension = '.tmp') {
    const tempPath = path.join(__dirname, `temp_${Date.now()}${extension}`);
    fs.writeFileSync(tempPath, content);
    this.tempFiles.push(tempPath);
    return tempPath;
  }

  cleanup() {
    this.tempFiles.forEach(file => {
      if (fs.existsSync(file)) {
        fs.unlinkSync(file);
      }
    });
    this.tempFiles = [];
  }
}

class TestAssertions {
  static expectValidJSON(jsonString) {
    let parsed;
    expect(() => {
      parsed = JSON.parse(jsonString);
    }).not.toThrow();
    return parsed;
  }

  static expectValidSMILES(smiles) {
    expect(typeof smiles).toBe('string');
    expect(smiles.length).toBeGreaterThan(0);
    // Basic SMILES validation - should contain chemical symbols
    expect(smiles).toMatch(/[CHNOPS]/);
  }

  static expectValidChemicalName(name) {
    expect(typeof name).toBe('string');
    expect(name.length).toBeGreaterThan(0);
    expect(name.trim()).toBe(name);
  }
}

// Test data helpers
const getTestMolecule = (type = 'simple') => {
  const molecules = {
    simple: {
      name: 'water',
      smiles: 'O',
      expectedAtoms: ['O', 'H', 'H']
    },
    complex: {
      name: 'caffeine', 
      smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
      expectedAtoms: ['C', 'N', 'O']
    },
    invalid: {
      name: 'invalid-molecule',
      smiles: 'INVALID',
      expectedAtoms: []
    }
  };
  return molecules[type] || molecules.simple;
};

const createTestRequest = (data = {}) => {
  return {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json'
    },
    body: JSON.stringify(data)
  };
};

module.exports = {
  TestFileManager,
  TestAssertions, 
  getTestMolecule,
  createTestRequest
}; 