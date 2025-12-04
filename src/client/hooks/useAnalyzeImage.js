/**
 * useAnalyzeImage - Predict molecular composition from images
 * 
 * This hook provides a simple interface for analyzing uploaded images
 * or photos to identify objects and their molecular composition.
 * 
 * @example
 * const { analyzeImage, isAnalyzing, molecules } = useAnalyzeImage();
 * 
 * // User uploads a photo of coffee
 * const handleFileUpload = async (file) => {
 *   const result = await analyzeImage(file);
 *   // Shows caffeine, chlorogenic acid, etc.
 * };
 */

import { useState, useCallback } from 'react';
import { useApi } from './useApi';  // Compose on top of existing hook

export function useAnalyzeImage() {
  const [isAnalyzing, setIsAnalyzing] = useState(false);
  const [molecules, setMolecules] = useState([]);
  const [identifiedObject, setIdentifiedObject] = useState('');
  const [error, setError] = useState(null);
  
  const { structuralizeImage } = useApi();  // Use existing functionality
  
  /**
   * Convert file to base64
   */
  const fileToBase64 = (file) => {
    return new Promise((resolve, reject) => {
      const reader = new FileReader();
      reader.onload = () => resolve(reader.result);
      reader.onerror = reject;
      reader.readAsDataURL(file);
    });
  };
  
  /**
   * Analyze an image file to find what's in it and its molecules
   * @param {File|string} image - Image file or base64 string
   * @param {Object} options - Analysis options
   * @param {number} options.x - X coordinate of interest point
   * @param {number} options.y - Y coordinate of interest point
   * @returns {Promise<{object: string, molecules: Array}>} What was found
   */
  const analyzeImage = useCallback(async (image, options = {}) => {
    setIsAnalyzing(true);
    setError(null);
    
    try {
      // Convert file to base64 if needed
      let imageData;
      if (image instanceof File) {
        imageData = await fileToBase64(image);
      } else if (typeof image === 'string') {
        imageData = image;
      } else {
        throw new Error('Invalid image format');
      }
      
      // Use existing API method with proper parameters
      const result = await structuralizeImage(
        imageData,
        options.x,
        options.y,
        options.x,  // cropMiddleX
        options.y,  // cropMiddleY
        options.cropSize || (options.x ? 300 : undefined)
      );
      
      // Transform to user-friendly format
      const molecules = (result.chemicals || result.molecules || []).map(chem => ({
        name: chem.name || chem.chemicalName,
        structure3D: chem.sdfPath,
        status: chem.status || 'ok'
      }));
      
      setMolecules(molecules);
      setIdentifiedObject(result.object || result.result?.object || 'Unknown object');
      
      return {
        object: result.object,
        molecules,
        boundingBox: result.recommendedBox
      };
      
    } catch (err) {
      const message = err.message || 'Failed to analyze image';
      setError(message);
      throw new Error(message);
    } finally {
      setIsAnalyzing(false);
    }
  }, [post]);
  
  /**
   * Analyze image from URL
   */
  const analyzeImageUrl = useCallback(async (url, options = {}) => {
    setIsAnalyzing(true);
    setError(null);
    
    try {
      // Use existing API method for URL
      const result = await structuralizeImage(
        url,
        options.x,
        options.y,
        options.x,
        options.y,
        options.cropSize
      );
      
      const molecules = (result.chemicals || result.molecules || []).map(chem => ({
        name: chem.name || chem.chemicalName,
        structure3D: chem.sdfPath,
        status: chem.status || 'ok'
      }));
      
      setMolecules(molecules);
      setIdentifiedObject(result.object || result.result?.object || 'Unknown object');
      
      return {
        object: result.object,
        molecules
      };
      
    } catch (err) {
      const message = err.message || 'Failed to analyze image from URL';
      setError(message);
      throw new Error(message);
    } finally {
      setIsAnalyzing(false);
    }
  }, [post]);
  
  /**
   * Clear results
   */
  const clear = useCallback(() => {
    setMolecules([]);
    setIdentifiedObject('');
    setError(null);
  }, []);
  
  return {
    // Actions
    analyzeImage,
    analyzeImageUrl,
    clear,
    
    // State
    isAnalyzing,
    molecules,
    identifiedObject,
    error,
    
    // Helpers
    hasMolecules: molecules.length > 0,
    moleculeCount: molecules.length
  };
}


