/**
 * useAnalyzeText - Predict molecular composition from text
 * 
 * This hook provides a simple interface for users who want to type something
 * and see what molecules it contains. It abstracts away all the complexity
 * of API calls, caching, and error handling.
 * 
 * @example
 * const { analyze, isAnalyzing, molecules, error } = useAnalyzeText();
 * 
 * // User types "coffee" and clicks analyze
 * const handleAnalyze = async () => {
 *   const result = await analyze('coffee');
 *   // result.molecules = [{name: 'Caffeine', formula: 'C8H10N4O2', sdfPath: '...'}]
 * };
 */

import { useState, useCallback } from 'react';
import { useApi } from './useApi';  // Compose on top of existing hook

export function useAnalyzeText() {
  const [isAnalyzing, setIsAnalyzing] = useState(false);
  const [molecules, setMolecules] = useState([]);
  const [error, setError] = useState(null);
  
  const { structuresFromText } = useApi();  // Use existing functionality
  
  /**
   * Analyze text to find molecular composition
   * @param {string} text - What to analyze (e.g., "coffee", "aspirin", "plastic bottle")
   * @param {Object} options - Analysis options
   * @returns {Promise<{molecules: Array, object: string}>} Analysis results
   */
  const analyze = useCallback(async (text, options = {}) => {
    setIsAnalyzing(true);
    setError(null);

    try {
      // Use existing API method
      const result = await structuresFromText(text.trim(), 'GPT-5');

      // Transform to user-friendly format
      const molecules = (result.chemicals || result.molecules || []).map(chem => ({
        name: chem.name || chem.chemicalName,
        formula: chem.formula || extractFormula(chem.name),
        structure3D: chem.sdfPath,
        status: chem.status || 'ok'
      }));

      setMolecules(molecules);

      return {
        molecules,
        object: result.object || text,
        source: result.source
      };

    } catch (err) {
      const message = err.message || 'Failed to analyze text';
      setError(message);
      throw new Error(message);
    } finally {
      setIsAnalyzing(false);
    }
  }, [post]);
  
  /**
   * Clear current results
   */
  const clear = useCallback(() => {
    setMolecules([]);
    setError(null);
  }, []);
  
  return {
    // Actions
    analyze,
    clear,
    
    // State
    isAnalyzing,
    molecules,
    error,
    
    // Helpers
    hasMolecules: molecules.length > 0,
    moleculeCount: molecules.length
  };
}

// Helper to extract formula (would be better from API)
function extractFormula(name) {
  // Common formulas - in real app, get from API
  const formulas = {
    'Caffeine': 'C₈H₁₀N₄O₂',
    'Water': 'H₂O',
    'Aspirin': 'C₉H₈O₄',
    'Glucose': 'C₆H₁₂O₆'
  };
  return formulas[name] || '';
}


