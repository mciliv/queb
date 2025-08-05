import React, { useState, useEffect, useRef, useCallback } from 'react';
import TextInput from './TextInput';
import ModeSelector from './ModeSelector';
import CameraSection from './CameraSection';
import PhotoSection from './PhotoSection';
import Results from './Results';
import PaymentSection from './PaymentSection';
import { usePayment } from './PaymentContext';
import { useApi } from '../hooks/useApi';

const MainLayout = ({ 
  isProcessing,
  setIsProcessing,
  viewers,
  setViewers,
  currentAnalysisType,
  setCurrentAnalysisType,
  lastAnalysis,
  setLastAnalysis
}) => {
  const [objectInput, setObjectInput] = useState('');
  const [cameraMode, setCameraMode] = useState(false);
  const [photoMode, setPhotoMode] = useState(false);
  const [showShortcuts, setShowShortcuts] = useState(false);
  const [error, setError] = useState('');
  const [retryCount, setRetryCount] = useState(0);
  const [lastSuccessfulAnalysis, setLastSuccessfulAnalysis] = useState(null);
  const { checkPaymentRequired } = usePayment();
  const { analyzeText, generateSDFs, error: apiError } = useApi();
  const maxRetries = 3;

  // Clear error when input changes
  useEffect(() => {
    if (error && objectInput) {
      setError('');
    }
  }, [objectInput, error]);

  // Handle API errors
  useEffect(() => {
    if (apiError) {
      setError(`API Error: ${apiError}`);
    }
  }, [apiError]);

  // Handle keyboard shortcuts
  useEffect(() => {
    const handleKeyDown = (event) => {
      // Ignore shortcuts when typing in input fields
      if (event.target.tagName === 'INPUT' || event.target.tagName === 'TEXTAREA') {
        return;
      }

      const isMac = navigator.platform.toUpperCase().indexOf('MAC') >= 0;
      const modifier = isMac ? event.metaKey : event.ctrlKey;

      switch (event.key.toLowerCase()) {
        case 'k':
          if (modifier) {
            event.preventDefault();
            document.getElementById('object-input')?.focus();
          }
          break;
        case 'c':
          if (modifier) {
            event.preventDefault();
            setCameraMode(prev => !prev);
            setPhotoMode(false);
          }
          break;
        case 'p':
          if (modifier) {
            event.preventDefault();
            setPhotoMode(prev => !prev);
            setCameraMode(false);
          }
          break;
        case 'escape':
          // Clear modes and focus
          setCameraMode(false);
          setPhotoMode(false);
          setShowShortcuts(false);
          setError('');
          document.getElementById('object-input')?.focus();
          break;
        case 'backspace':
          if (modifier) {
            event.preventDefault();
            setViewers([]);
            setLastAnalysis(null);
            setError('');
          }
          break;
        case 'enter':
          if (modifier && objectInput.trim()) {
            event.preventDefault();
            handleTextAnalysis(objectInput);
          }
          break;
        case '?':
          if (modifier) {
            event.preventDefault();
            setShowShortcuts(prev => !prev);
          }
          break;
      }
    };

    document.addEventListener('keydown', handleKeyDown);
    return () => document.removeEventListener('keydown', handleKeyDown);
  }, [objectInput, setViewers, setLastAnalysis]);

  const handleTextAnalysis = useCallback(async (value) => {
    if (isProcessing || !value.trim()) return;

    if (checkPaymentRequired()) {
      // Show payment modal
      return;
    }

    setIsProcessing(true);
    setCurrentAnalysisType('text');
    setError('');
    
    try {
      const result = await analyzeText(value);
      
      // Handle both 'molecules' and 'chemicals' field names for compatibility
      const molecules = result.molecules || result.chemicals || [];
      
      if (molecules && molecules.length > 0) {
        // Extract SMILES strings for SDF generation
        const smilesArray = molecules.map(mol => mol.smiles).filter(Boolean);
        
        if (smilesArray.length > 0) {
          try {
            // Generate SDF data from SMILES
            const sdfResult = await generateSDFs(smilesArray, false);
            
            // Create viewers with SDF data
            const newViewers = molecules.map((mol, index) => {
              const sdfPath = sdfResult.sdfPaths && sdfResult.sdfPaths[index];
              return {
                name: mol.name || value,
                sdfData: sdfPath ? `file://${sdfPath}` : null, // Use file path for SDF data
                smiles: mol.smiles
              };
            });
            
            setViewers(prev => [...prev, ...newViewers]);
            setLastSuccessfulAnalysis(result);
            setRetryCount(0); // Reset retry count on success
          } catch (sdfError) {
            console.error('SDF generation failed:', sdfError);
            // Still show molecules even if SDF generation fails
            const newViewers = molecules.map(mol => ({
              name: mol.name || value,
              sdfData: null,
              smiles: mol.smiles
            }));
            setViewers(prev => [...prev, ...newViewers]);
            setLastSuccessfulAnalysis(result);
            setRetryCount(0);
          }
        } else {
          setError('No valid SMILES found in the analysis results.');
        }
      } else {
        setError('No molecules found for this input. Try a different chemical name or formula.');
      }
      
      setLastAnalysis(result);
      setObjectInput('');
    } catch (error) {
      console.error('Analysis failed:', error);
      
      // Implement retry logic
      if (retryCount < maxRetries) {
        setRetryCount(prev => prev + 1);
        setError(`Analysis failed. Retrying... (${retryCount + 1}/${maxRetries})`);
        
        // Retry after a short delay
        setTimeout(() => {
          handleTextAnalysis(value);
        }, 1000 * (retryCount + 1)); // Exponential backoff
      } else {
        setError(`Analysis failed after ${maxRetries} attempts. Please check your input and try again.`);
        setRetryCount(0);
      }
    } finally {
      setIsProcessing(false);
    }
  }, [isProcessing, analyzeText, generateSDFs, setViewers, setLastAnalysis, checkPaymentRequired, retryCount]);

  const handleAnalysisComplete = useCallback(async (result) => {
    // Handle both 'molecules' and 'chemicals' field names for compatibility
    const molecules = result.molecules || result.chemicals || [];
    
    if (molecules && molecules.length > 0) {
      // Extract SMILES strings for SDF generation
      const smilesArray = molecules.map(mol => mol.smiles).filter(Boolean);
      
      if (smilesArray.length > 0) {
        try {
          // Generate SDF data from SMILES
          const sdfResult = await generateSDFs(smilesArray, false);
          
          // Create viewers with SDF data
          const newViewers = molecules.map((mol, index) => {
            const sdfPath = sdfResult.sdfPaths && sdfResult.sdfPaths[index];
            return {
              name: mol.name || 'Captured object',
              sdfData: sdfPath ? `file://${sdfPath}` : null,
              smiles: mol.smiles
            };
          });
          
          setViewers(prev => [...prev, ...newViewers]);
          setLastSuccessfulAnalysis(result);
        } catch (sdfError) {
          console.error('SDF generation failed:', sdfError);
          // Still show molecules even if SDF generation fails
          const newViewers = molecules.map(mol => ({
            name: mol.name || 'Captured object',
            sdfData: null,
            smiles: mol.smiles
          }));
          setViewers(prev => [...prev, ...newViewers]);
          setLastSuccessfulAnalysis(result);
        }
      } else {
        // Show molecules without SDF data
        const newViewers = molecules.map(mol => ({
          name: mol.name || 'Captured object',
          sdfData: null,
          smiles: mol.smiles
        }));
        setViewers(prev => [...prev, ...newViewers]);
        setLastSuccessfulAnalysis(result);
      }
    }
    
    setLastAnalysis(result);
    setError('');
  }, [setViewers, setLastAnalysis, generateSDFs]);