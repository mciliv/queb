/**
 * Unit tests for visual test helper utilities
 */

const {
  sanitizeSmiles,
  normalizePath,
  resolveSdfPath,
  createViewerObject,
  processSmilesToViewers
} = require('../../../src/client/utils/visual-test-helpers.js');

describe('Visual Test Helpers', () => {
  
  describe('sanitizeSmiles', () => {
    it('should replace special characters with underscores', () => {
      expect(sanitizeSmiles('C(=O)O')).toBe('C___O_O');
    });

    it('should replace equals sign with double underscores', () => {
      expect(sanitizeSmiles('C=O')).toBe('C__O');
    });

    it('should handle complex SMILES notation', () => {
      expect(sanitizeSmiles('CN1C=NC2=C1C(=O)N(C(=O)N2C)C')).toBe(
        'CN1C__NC2__C1C___O_N_C___O_N2C_C'
      );
    });

    it('should preserve alphanumeric characters', () => {
      expect(sanitizeSmiles('CCO123')).toBe('CCO123');
    });

    it('should handle empty string', () => {
      expect(sanitizeSmiles('')).toBe('');
    });
  });

  describe('normalizePath', () => {
    it('should add leading slash to path without one', () => {
      expect(normalizePath('sdf_files/test.sdf')).toBe('/sdf_files/test.sdf');
    });

    it('should normalize multiple leading slashes to one', () => {
      expect(normalizePath('///sdf_files/test.sdf')).toBe('/sdf_files/test.sdf');
    });

    it('should keep single leading slash unchanged', () => {
      expect(normalizePath('/sdf_files/test.sdf')).toBe('/sdf_files/test.sdf');
    });

    it('should handle empty string', () => {
      expect(normalizePath('')).toBe('');
    });
  });

  describe('resolveSdfPath', () => {
    it('should return deterministic path when found in returned paths', () => {
      const smiles = 'CCO';
      const returnedPaths = ['/sdf_files/CCO.sdf', '/sdf_files/other.sdf'];
      
      expect(resolveSdfPath(smiles, 0, returnedPaths)).toBe('/sdf_files/CCO.sdf');
    });

    it('should return index-based path when deterministic not found', () => {
      const smiles = 'CCO';
      const returnedPaths = ['/sdf_files/random123.sdf', '/sdf_files/other.sdf'];
      
      expect(resolveSdfPath(smiles, 0, returnedPaths)).toBe('/sdf_files/random123.sdf');
    });

    it('should fallback to deterministic path when no matches', () => {
      const smiles = 'C=O';
      const returnedPaths = ['/sdf_files/other.sdf'];
      
      expect(resolveSdfPath(smiles, 5, returnedPaths)).toBe('/sdf_files/C__O.sdf');
    });

    it('should handle empty returnedPaths', () => {
      const smiles = 'CCO';
      
      expect(resolveSdfPath(smiles, 0, [])).toBe('/sdf_files/CCO.sdf');
    });

    it('should normalize paths before comparison', () => {
      const smiles = 'CCO';
      const returnedPaths = ['//sdf_files/CCO.sdf']; // Double slash
      
      // Should match despite different slash count
      expect(resolveSdfPath(smiles, 0, returnedPaths)).toBe('/sdf_files/CCO.sdf');
    });
  });

  describe('createViewerObject', () => {
    it('should create viewer object with name from map', () => {
      const nameMap = { 'CCO': 'Ethanol' };
      const result = createViewerObject('CCO', '/sdf_files/CCO.sdf', nameMap);
      
      expect(result).toEqual({
        name: 'Ethanol',
        sdfData: 'file:///sdf_files/CCO.sdf',
        smiles: 'CCO'
      });
    });

    it('should use SMILES as name when not in map', () => {
      const result = createViewerObject('CCO', '/sdf_files/CCO.sdf', {});
      
      expect(result).toEqual({
        name: 'CCO',
        sdfData: 'file:///sdf_files/CCO.sdf',
        smiles: 'CCO'
      });
    });

    it('should handle null sdfPath', () => {
      const result = createViewerObject('CCO', null, {});
      
      expect(result).toEqual({
        name: 'CCO',
        sdfData: null,
        smiles: 'CCO'
      });
    });

    it('should work with empty nameMap', () => {
      const result = createViewerObject('O', '/sdf_files/O.sdf');
      
      expect(result.name).toBe('O');
      expect(result.smiles).toBe('O');
    });
  });

  describe('processSmilesToViewers', () => {
    const nameMap = {
      'O': 'Water',
      'CCO': 'Ethanol',
      'C=O': 'Formaldehyde'
    };

    it('should process array of SMILES into viewer objects', () => {
      const smilesArray = ['O', 'CCO'];
      const returnedPaths = ['/sdf_files/O.sdf', '/sdf_files/CCO.sdf'];
      
      const result = processSmilesToViewers(smilesArray, returnedPaths, nameMap);
      
      expect(result).toHaveLength(2);
      expect(result[0]).toEqual({
        name: 'Water',
        sdfData: 'file:///sdf_files/O.sdf',
        smiles: 'O'
      });
      expect(result[1]).toEqual({
        name: 'Ethanol',
        sdfData: 'file:///sdf_files/CCO.sdf',
        smiles: 'CCO'
      });
    });

    it('should handle empty arrays', () => {
      const result = processSmilesToViewers([], [], {});
      expect(result).toEqual([]);
    });

    it('should handle mismatched array lengths', () => {
      const smilesArray = ['O', 'CCO', 'C=O'];
      const returnedPaths = ['/sdf_files/O.sdf']; // Only one path
      
      const result = processSmilesToViewers(smilesArray, returnedPaths, nameMap);
      
      expect(result).toHaveLength(3);
      expect(result[0].sdfData).toBe('file:///sdf_files/O.sdf');
      // Others should fallback to deterministic paths
      expect(result[1].sdfData).toBe('file:///sdf_files/CCO.sdf');
      expect(result[2].sdfData).toBe('file:///sdf_files/C__O.sdf');
    });

    it('should work without nameMap or returnedPaths', () => {
      const smilesArray = ['O', 'CCO'];
      
      const result = processSmilesToViewers(smilesArray);
      
      expect(result).toHaveLength(2);
      expect(result[0].name).toBe('O'); // Using SMILES as name
      expect(result[0].sdfData).toBe('file:///sdf_files/O.sdf');
    });

    it('should handle complex SMILES with special characters', () => {
      const smilesArray = ['CN1C=NC2=C1C(=O)N(C(=O)N2C)C'];
      const result = processSmilesToViewers(smilesArray, [], {});
      
      expect(result[0].sdfData).toContain('CN1C__NC2__C1C___O_N_C___O_N2C_C');
    });
  });

  describe('Integration scenarios', () => {
    it('should handle a complete visual test scenario', () => {
      // Simulating a "Coffee" test with caffeine, water, ethanol
      const test = {
        label: 'Coffee',
        smilesList: ['CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'O', 'CCO']
      };
      
      const serverResponse = {
        sdfPaths: [
          '/sdf_files/CN1C__NC2__C1C___O_N_C___O_N2C_C.sdf',
          '/sdf_files/O.sdf',
          '/sdf_files/CCO.sdf'
        ]
      };
      
      const nameMap = {
        'CN1C=NC2=C1C(=O)N(C(=O)N2C)C': 'Caffeine',
        'O': 'Water',
        'CCO': 'Ethanol'
      };
      
      const viewers = processSmilesToViewers(
        test.smilesList,
        serverResponse.sdfPaths,
        nameMap
      );
      
      expect(viewers).toHaveLength(3);
      expect(viewers[0].name).toBe('Caffeine');
      expect(viewers[1].name).toBe('Water');
      expect(viewers[2].name).toBe('Ethanol');
      
      viewers.forEach(viewer => {
        expect(viewer.sdfData).toMatch(/^file:\/\//);
        expect(viewer.smiles).toBeTruthy();
      });
    });

    it('should gracefully handle server errors with fallback paths', () => {
      const smilesArray = ['O', 'CCO'];
      const emptyResponse = []; // Server returned no paths
      const nameMap = { 'O': 'Water', 'CCO': 'Ethanol' };
      
      const viewers = processSmilesToViewers(smilesArray, emptyResponse, nameMap);
      
      expect(viewers).toHaveLength(2);
      // Should still have valid fallback paths
      expect(viewers[0].sdfData).toBe('file:///sdf_files/O.sdf');
      expect(viewers[1].sdfData).toBe('file:///sdf_files/CCO.sdf');
      // Names should still be from map
      expect(viewers[0].name).toBe('Water');
      expect(viewers[1].name).toBe('Ethanol');
    });
  });
});

