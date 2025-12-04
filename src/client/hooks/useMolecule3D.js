/**
 * useMolecule3D - Generate and manage 3D molecular structures
 * 
 * This hook focuses on the 3D visualization aspect - taking molecular
 * data and making it viewable in 3D. It handles structure generation
 * and provides data ready for 3D viewers.
 * 
 * @example
 * const { generate3D, structures, isGenerating } = useMolecule3D();
 * 
 * // After analysis finds caffeine
 * const structure = await generate3D('caffeine');
 * // Returns path to 3D structure file ready for viewing
 */

import { useState, useCallback } from 'react';
import { useApi } from './useApi';  // Compose on top of existing hook

export function useMolecule3D() {
  const [isGenerating, setIsGenerating] = useState(false);
  const [structures, setStructures] = useState(new Map());
  const [error, setError] = useState(null);
  
  const { nameToSdf, generateSDFs } = useApi();  // Use existing functionality
  
  /**
   * Generate 3D structure from chemical name
   * @param {string} chemicalName - Name of chemical (e.g., "caffeine")
   * @param {Object} options - Generation options
   * @param {boolean} options.forceRegenerate - Regenerate even if cached
   * @returns {Promise<{name: string, path: string, format: string}>} 3D structure info
   */
  const generate3DFromName = useCallback(async (chemicalName, options = {}) => {
    if (!chemicalName) {
      throw new Error('Chemical name is required');
    }
    
    // Check cache first
    if (!options.forceRegenerate && structures.has(chemicalName)) {
      return structures.get(chemicalName);
    }
    
    setIsGenerating(true);
    setError(null);
    
    try {
      const result = await nameToSdf(
        chemicalName,
        options.forceRegenerate || false
      );
      
      const structure = {
        name: chemicalName,
        path: result.sdfPath,
        format: 'sdf',
        smiles: result.smiles
      };
      
      // Cache the result
      setStructures(prev => new Map(prev).set(chemicalName, structure));
      
      return structure;
      
    } catch (err) {
      const message = `Failed to generate 3D structure for ${chemicalName}`;
      setError(message);
      throw new Error(message);
    } finally {
      setIsGenerating(false);
    }
  }, [structures, post]);
  
  /**
   * Generate 3D structures for multiple molecules
   * @param {Array<{name: string, smiles?: string}>} molecules - Molecules to generate
   * @returns {Promise<Array>} Generated structures
   */
  const generateMultiple3D = useCallback(async (molecules) => {
    if (!molecules || molecules.length === 0) {
      return [];
    }
    
    setIsGenerating(true);
    setError(null);
    
    try {
      // Extract SMILES if available, otherwise we'll use names
      const smilesArray = molecules
        .filter(mol => mol.smiles)
        .map(mol => mol.smiles);
      
      let paths = [];
      
      if (smilesArray.length > 0) {
        // Generate from SMILES (more accurate)
        const result = await generateSDFs(smilesArray, false);
        paths = result.sdfPaths || [];
      }
      
      // For molecules without SMILES, generate by name
      const withoutSmiles = molecules.filter(mol => !mol.smiles);
      for (const mol of withoutSmiles) {
        try {
          const result = await generate3DFromName(mol.name);
          paths.push(result.path);
        } catch (err) {
          console.warn(`Failed to generate structure for ${mol.name}`);
        }
      }
      
      // Create structure objects
      const structures = molecules.map((mol, index) => ({
        name: mol.name,
        path: paths[index] || null,
        format: 'sdf',
        smiles: mol.smiles
      }));
      
      return structures.filter(s => s.path); // Only return successful ones
      
    } catch (err) {
      const message = 'Failed to generate some 3D structures';
      setError(message);
      throw new Error(message);
    } finally {
      setIsGenerating(false);
    }
  }, [post, generate3DFromName]);
  
  /**
   * Get viewer-ready data for a molecule
   * @param {string} structurePath - Path to structure file
   * @returns {Object} Data ready for 3D viewer component
   */
  const getViewerData = useCallback((structurePath) => {
    if (!structurePath) return null;
    
    return {
      sdfData: `file://${structurePath}`,
      format: 'sdf',
      viewerOptions: {
        backgroundColor: '#000000',
        style: 'ball-and-stick',
        colorScheme: 'rasmol'
      }
    };
  }, []);
  
  /**
   * Clear all cached structures
   */
  const clearCache = useCallback(() => {
    setStructures(new Map());
    setError(null);
  }, []);
  
  return {
    // Actions
    generate3DFromName,
    generateMultiple3D,
    getViewerData,
    clearCache,
    
    // State
    isGenerating,
    structures: Array.from(structures.values()),
    error,
    
    // Helpers
    hasStructures: structures.size > 0,
    getStructure: (name) => structures.get(name)
  };
}


