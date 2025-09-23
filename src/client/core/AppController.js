/**
 * Application Controller
 * 
 * This module applies the "Deep Modules" principle by separating business logic
 * from UI concerns, providing a clean interface for application state management.
 * 
 * Philosophy: "Classes should be deep: they should have simple interfaces that hide complex implementations"
 */

import { useState, useCallback, useRef } from 'react';
import { useApi } from '../hooks/useApi';
import { SMILES_NAME_MAP } from './constants.js';
import logger from './logger.js';

/**
 * Custom hook for application state management
 * Deep module: complex state logic hidden behind simple interface
 */
export function useAppController() {
  // Core state
  const [columns, setColumns] = useState([]);
  const [isProcessing, setIsProcessing] = useState(false);
  const [error, setError] = useState('');
  
  // Input modes
  const [activeMode, setActiveMode] = useState(null); // 'camera', 'photo', 'link', or null
  const [objectInput, setObjectInput] = useState('');
  
  // Settings
  const [columnMode, setColumnMode] = useState('replace'); // 'replace' or 'accumulate'
  const [showSettings, setShowSettings] = useState(false);
  
  // API integration
  const { structuresFromText: analyzeText, generateSDFs, twoNamesToSdf } = useApi();
  
  // Processing state tracking
  const processingRef = useRef(false);

  /**
   * INPUT HANDLING - Simple interface for complex input processing
   */

  const handleTextInput = useCallback(async (text) => {
    if (!text?.trim() || processingRef.current) return;
    
    processingRef.current = true;
    setIsProcessing(true);
    setError('');
    
    try {
      await processTextAnalysis(text.trim());
      setObjectInput(''); // Clear input on success
    } catch (error) {
      handleError(error, 'text analysis');
    } finally {
      processingRef.current = false;
      setIsProcessing(false);
    }
  }, []);

  const handleImageAnalysis = useCallback(async (result) => {
    if (processingRef.current) return;
    
    processingRef.current = true;
    setIsProcessing(true);
    setError('');
    
    try {
      await processImageAnalysis(result);
    } catch (error) {
      handleError(error, 'image analysis');
    } finally {
      processingRef.current = false;
      setIsProcessing(false);
    }
  }, []);

  /**
   * MODE MANAGEMENT - Simple interface for input mode switching
   */

  const switchMode = useCallback((mode) => {
    // Clear other modes when switching
    setActiveMode(current => current === mode ? null : mode);
    setError(''); // Clear any existing errors
  }, []);

  const clearActiveMode = useCallback(() => {
    setActiveMode(null);
  }, []);

  /**
   * COLUMN MANAGEMENT - Simple interface for result management
   */

  const removeColumn = useCallback((columnId) => {
    setColumns(prev => prev.filter(col => col.id !== columnId));
  }, []);

  const clearAllColumns = useCallback(() => {
    setColumns([]);
  }, []);

  /**
   * SETTINGS MANAGEMENT
   */

  const toggleSettings = useCallback(() => {
    setShowSettings(prev => !prev);
  }, []);

  const updateColumnMode = useCallback((mode) => {
    setColumnMode(mode);
    logger.info('Column mode changed', { mode });
  }, []);

  /**
   * PRIVATE IMPLEMENTATION - Complex logic hidden from components
   */

  /**
   * Process text analysis with intelligent routing
   */
  const processTextAnalysis = async (text) => {
    logger.info(`🔍 Starting text analysis for: "${text}"`);
    
    const columnId = Date.now();
    const newColumn = createColumn(columnId, text);
    
    // Add column based on mode
    if (columnMode === 'replace') {
      setColumns([newColumn]);
    } else {
      setColumns(prev => [...prev, newColumn]);
    }
    
    try {
      // Intelligent analysis routing
      const result = await routeTextAnalysis(text);
      const processedMolecules = await processMolecules(result);
      
      // Update column with results
      updateColumn(columnId, {
        viewers: processedMolecules,
        loading: false,
        failed: false
      });
      
      logger.info(`✅ Text analysis complete: Found ${processedMolecules.length} molecules`);
      
    } catch (error) {
      updateColumn(columnId, { loading: false, failed: true });
      throw error;
    }
  };

  /**
   * Process image analysis results
   */
  const processImageAnalysis = async (analysisResult) => {
    const molecules = analysisResult?.molecules || analysisResult?.chemicals || [];
    const objectLabel = getObjectLabel(analysisResult);
    
    const columnId = Date.now();
    const newColumn = createColumn(columnId, objectLabel);
    
    if (columnMode === 'replace') {
      setColumns([newColumn]);
    } else {
      setColumns(prev => [...prev, newColumn]);
    }
    
    try {
      if (molecules.length > 0) {
        const processedMolecules = await processMolecules({ molecules });
        updateColumn(columnId, {
          viewers: processedMolecules,
          loading: false,
          failed: false
        });
      } else {
        updateColumn(columnId, {
          viewers: [],
          loading: false,
          failed: true
        });
      }
    } catch (error) {
      updateColumn(columnId, { loading: false, failed: true });
      throw error;
    }
  };

  /**
   * Route text analysis to appropriate service
   */
  const routeTextAnalysis = async (text) => {
    // Check if it's exactly two comma-separated names
    const parts = text.split(',').map(s => s.trim()).filter(Boolean);
    
    if (parts.length === 2) {
      try {
        const response = await twoNamesToSdf(parts, false);
        return { object: text, molecules: response.molecules || [] };
      } catch (error) {
        // Fallback to standard analysis
        logger.warn('Two-names analysis failed, falling back to standard analysis');
      }
    }
    
    // Standard text analysis
    return await analyzeText(text);
  };

  /**
   * Process molecules and generate 3D structures
   */
  const processMolecules = async (result) => {
    const molecules = result.molecules || result.chemicals || [];
    
    // Check for pre-computed SDFs first
    const precomputed = molecules.filter(m => m.sdfPath);
    if (precomputed.length > 0) {
      return precomputed.map(m => ({
        name: m.name || result.object,
        sdfData: m.sdfPath.startsWith('file://') ? m.sdfPath : `file://${m.sdfPath}`,
        smiles: m.smiles
      }));
    }
    
    // Generate SDFs from SMILES
    const smilesArray = molecules.map(mol => mol.smiles).filter(Boolean);
    if (smilesArray.length === 0) {
      return [];
    }
    
    try {
      const sdfResult = await generateSDFs(smilesArray, false);
      const smilesToSdf = new Map();
      sdfResult.sdfPaths?.forEach((path, index) => {
        smilesToSdf.set(smilesArray[index], path);
      });
      
      return molecules.map(mol => {
        const sdfPath = smilesToSdf.get(mol.smiles);
        return {
          name: mol.name || SMILES_NAME_MAP[mol.smiles] || mol.smiles || result.object,
          sdfData: sdfPath ? `file://${sdfPath}` : null,
          smiles: mol.smiles
        };
      });
      
    } catch (error) {
      logger.error('SDF generation failed:', error);
      throw new Error('Failed to generate molecular structures');
    }
  };

  /**
   * Utility functions
   */
  const createColumn = (id, query) => ({
    id,
    query,
    viewers: [],
    loading: true,
    failed: false
  });

  const updateColumn = (id, updates) => {
    setColumns(prev => prev.map(col => 
      col.id === id ? { ...col, ...updates } : col
    ));
  };

  const getObjectLabel = (result) => {
    return (result && (result.object || (result.result && result.result.object))) || 
           (activeMode === 'camera' ? 'Camera capture' : 'Image capture');
  };

  const handleError = (error, context) => {
    const message = error?.message || 'An unknown error occurred';
    logger.error(`${context} failed:`, error);
    setError(`${context} failed: ${message}`);
  };

  /**
   * PUBLIC INTERFACE - Simple methods for components
   */
  return {
    // State
    columns,
    isProcessing,
    error,
    objectInput,
    activeMode,
    columnMode,
    showSettings,
    
    // Input handlers
    handleTextInput,
    handleImageAnalysis,
    setObjectInput,
    
    // Mode management
    switchMode,
    clearActiveMode,
    
    // Column management
    removeColumn,
    clearAllColumns,
    
    // Settings
    toggleSettings,
    updateColumnMode,
    
    // Utilities
    clearError: () => setError(''),
    
    // Mode checkers
    isCameraMode: activeMode === 'camera',
    isPhotoMode: activeMode === 'photo',
    isLinkMode: activeMode === 'link'
  };
}


