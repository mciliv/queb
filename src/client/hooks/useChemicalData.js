/**
 * useChemicalData - Look up detailed information about chemicals
 * 
 * This hook provides a focused interface for getting chemical properties,
 * formulas, names, and other data. It's about learning more about a
 * specific chemical rather than analyzing objects.
 * 
 * @example
 * const { lookupChemical, chemicalInfo } = useChemicalData();
 * 
 * // User wants to know more about caffeine
 * const info = await lookupChemical('caffeine');
 * // Returns formula, IUPAC name, properties, etc.
 */

import { useState, useCallback } from 'react';
import { useApi } from './useApi';  // Compose on top of existing hook

export function useChemicalData() {
  const [isLoading, setIsLoading] = useState(false);
  const [chemicalInfo, setChemicalInfo] = useState(null);
  const [error, setError] = useState(null);
  
  const { nameToSmiles } = useApi();  // Use existing functionality
  
  /**
   * Look up detailed chemical information
   * @param {string} chemicalName - Name of chemical to look up
   * @returns {Promise<Object>} Chemical properties and data
   */
  const lookupChemical = useCallback(async (chemicalName) => {
    if (!chemicalName?.trim()) {
      throw new Error('Chemical name is required');
    }
    
    setIsLoading(true);
    setError(null);
    
    try {
      // First resolve the name to get structured data
      const resolved = await nameToSmiles(chemicalName.trim());
      
      // Transform to user-friendly format
      const info = {
        // Basic identification
        name: resolved.title || chemicalName,
        iupacName: resolved.iupac || null,
        commonNames: resolved.synonyms || [],
        
        // Structure data
        smiles: resolved.smiles,
        formula: extractMolecularFormula(resolved.smiles),
        pubchemId: resolved.cid || null,
        
        // Properties (would be extended with more API calls)
        properties: {
          molecularWeight: calculateMolecularWeight(resolved.smiles),
          appearance: null, // Would come from extended data
          meltingPoint: null,
          boilingPoint: null
        },
        
        // Links for more info
        references: {
          pubchem: resolved.cid ? `https://pubchem.ncbi.nlm.nih.gov/compound/${resolved.cid}` : null,
          wikipedia: null // Could add Wikipedia API lookup
        }
      };
      
      setChemicalInfo(info);
      return info;
      
    } catch (err) {
      const message = `Failed to look up ${chemicalName}`;
      setError(message);
      throw new Error(message);
    } finally {
      setIsLoading(false);
    }
  }, [post]);
  
  /**
   * Look up chemical by CAS number
   */
  const lookupByCAS = useCallback(async (casNumber) => {
    // CAS format: XXX-XX-X
    if (!casNumber?.match(/^\d{2,7}-\d{2}-\d$/)) {
      throw new Error('Invalid CAS number format');
    }
    
    // Would implement CAS lookup
    return lookupChemical(casNumber);
  }, [lookupChemical]);
  
  /**
   * Look up chemical by SMILES
   */
  const lookupBySMILES = useCallback(async (smiles) => {
    if (!smiles?.trim()) {
      throw new Error('SMILES notation is required');
    }
    
    setIsLoading(true);
    setError(null);
    
    try {
      // For SMILES, we already have the structure
      const info = {
        name: 'Unknown compound',
        smiles: smiles,
        formula: extractMolecularFormula(smiles),
        properties: {
          molecularWeight: calculateMolecularWeight(smiles)
        }
      };
      
      setChemicalInfo(info);
      return info;
      
    } catch (err) {
      const message = 'Failed to process SMILES notation';
      setError(message);
      throw new Error(message);
    } finally {
      setIsLoading(false);
    }
  }, []);
  
  /**
   * Clear current data
   */
  const clear = useCallback(() => {
    setChemicalInfo(null);
    setError(null);
  }, []);
  
  return {
    // Actions
    lookupChemical,
    lookupByCAS,
    lookupBySMILES,
    clear,
    
    // State
    isLoading,
    chemicalInfo,
    error,
    
    // Helpers
    hasData: chemicalInfo !== null
  };
}

// Helper functions (simplified - real implementation would use cheminformatics library)
function extractMolecularFormula(smiles) {
  // Very simplified - real implementation would parse SMILES properly
  const formulas = {
    'CN1C=NC2=C1C(=O)N(C(=O)N2C)C': 'C₈H₁₀N₄O₂', // Caffeine
    'CC(=O)OC1=CC=CC=C1C(=O)O': 'C₉H₈O₄', // Aspirin
    'O': 'H₂O', // Water
    'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O': 'C₆H₁₂O₆' // Glucose
  };
  return formulas[smiles] || 'Unknown';
}

function calculateMolecularWeight(smiles) {
  // Very simplified - real implementation would calculate from SMILES
  const weights = {
    'CN1C=NC2=C1C(=O)N(C(=O)N2C)C': 194.19, // Caffeine
    'CC(=O)OC1=CC=CC=C1C(=O)O': 180.16, // Aspirin
    'O': 18.015, // Water
    'C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O': 180.16 // Glucose
  };
  return weights[smiles] || null;
}


