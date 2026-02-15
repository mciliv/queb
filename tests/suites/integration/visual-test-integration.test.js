/**
 * Integration tests for visual test system
 * Tests the complete flow of SDF generation and viewer object creation
 */

const { runVisualTest, runAllVisualTests } = require('../visual/visual-test-runner.js');
const { processSmilesToViewers } = require('../../../src/client/utils/visual-test-helpers.js');

describe('Visual Test Integration', () => {
  
  describe('runVisualTest', () => {
    it('should successfully run a valid test', async () => {
      const test = {
        label: 'Test Water',
        smilesList: ['O']
      };
      
      const result = await runVisualTest(test);
      
      expect(result.success).toBe(true);
      expect(result.viewers).toHaveLength(1);
      expect(result.viewers[0].smiles).toBe('O');
      expect(result.viewers[0].sdfData).toContain('file://');
    });

    it('should handle test with multiple molecules', async () => {
      const test = {
        label: 'Test Multiple',
        smilesList: ['O', 'CCO', 'C=O']
      };
      
      const result = await runVisualTest(test);
      
      expect(result.success).toBe(true);
      expect(result.viewers).toHaveLength(3);
      result.viewers.forEach(viewer => {
        expect(viewer.smiles).toBeTruthy();
        expect(viewer.name).toBeTruthy();
        expect(viewer.sdfData).toMatch(/^file:\/\//);
      });
    });

    it('should fail gracefully for test with no SMILES', async () => {
      const test = {
        label: 'Empty Test',
        smilesList: []
      };
      
      const result = await runVisualTest(test);
      
      expect(result.success).toBe(false);
      expect(result.reason).toBe('no_smiles');
    });

    it('should handle missing smilesList property', async () => {
      const test = {
        label: 'Invalid Test'
      };
      
      const result = await runVisualTest(test);
      
      expect(result.success).toBe(false);
    });
  });

  describe('runAllVisualTests', () => {
    it('should run all preset tests and return summary', async () => {
      const summary = await runAllVisualTests();
      
      expect(summary).toHaveProperty('passed');
      expect(summary).toHaveProperty('failed');
      expect(summary).toHaveProperty('results');
      expect(Array.isArray(summary.results)).toBe(true);
      expect(summary.results.length).toBeGreaterThan(0);
    });

    it('should have all tests pass with valid data', async () => {
      const summary = await runAllVisualTests();
      
      // All preset tests should be valid
      expect(summary.passed).toBeGreaterThan(0);
      expect(summary.failed).toBe(0);
    });
  });

  describe('End-to-end viewer generation', () => {
    it('should generate valid viewers for Coffee test', async () => {
      const mockServerResponse = {
        sdfPaths: [
          '/sdf_files/CN1C__NC2__C1C___O_N_C___O_N2C_C.sdf',
          '/sdf_files/O.sdf',
          '/sdf_files/CCO.sdf'
        ]
      };
      
      const smilesList = ['CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'O', 'CCO'];
      const nameMap = {
        'CN1C=NC2=C1C(=O)N(C(=O)N2C)C': 'Caffeine',
        'O': 'Water',
        'CCO': 'Ethanol'
      };
      
      const viewers = processSmilesToViewers(
        smilesList,
        mockServerResponse.sdfPaths,
        nameMap
      );
      
      expect(viewers).toHaveLength(3);
      expect(viewers[0].name).toBe('Caffeine');
      expect(viewers[1].name).toBe('Water');
      expect(viewers[2].name).toBe('Ethanol');
      
      // All should have valid file paths
      viewers.forEach(viewer => {
        expect(viewer.sdfData).toMatch(/^file:\/\//);
        expect(viewer.smiles).toBeTruthy();
      });
    });

    it('should handle partial server response gracefully', async () => {
      const partialResponse = {
        sdfPaths: ['/sdf_files/O.sdf'] // Only first one returned
      };
      
      const smilesList = ['O', 'CCO', 'C=O'];
      
      const viewers = processSmilesToViewers(
        smilesList,
        partialResponse.sdfPaths,
        {}
      );
      
      expect(viewers).toHaveLength(3);
      // All should still have valid paths (fallback to deterministic)
      viewers.forEach(viewer => {
        expect(viewer.sdfData).toBeTruthy();
        expect(viewer.sdfData).toMatch(/^file:\/\//);
      });
    });
  });
});

