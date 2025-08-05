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
  const { analyzeText, error: apiError } = useApi();
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
      
      if (result.molecules && result.molecules.length > 0) {
        const newViewers = result.molecules.map(mol => ({
          name: mol.name || value,
          sdfData: mol.sdf_data,
          smiles: mol.smiles
        }));
        setViewers(prev => [...prev, ...newViewers]);
        setLastSuccessfulAnalysis(result);
        setRetryCount(0); // Reset retry count on success
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
  }, [isProcessing, analyzeText, setViewers, setLastAnalysis, checkPaymentRequired, retryCount]);

  const handleAnalysisComplete = useCallback((result) => {
    if (result.molecules && result.molecules.length > 0) {
      const newViewers = result.molecules.map(mol => ({
        name: mol.name || 'Captured object',
        sdfData: mol.sdf_data,
        smiles: mol.smiles
      }));
      setViewers(prev => [...prev, ...newViewers]);
      setLastSuccessfulAnalysis(result);
    }
    
    setLastAnalysis(result);
    setError('');
  }, [setViewers, setLastAnalysis]);

  const handleRetry = useCallback(() => {
    if (lastSuccessfulAnalysis) {
      handleTextAnalysis(objectInput || lastSuccessfulAnalysis.query);
    }
  }, [lastSuccessfulAnalysis, objectInput, handleTextAnalysis]);

  return (
    <div className="app-container">
      <div className="main-app-interface">
        <div className="main-content-layout">
          {/* Left side: Analysis section */}
          <div className="analysis-section">
            <TextInput 
              value={objectInput}
              onChange={setObjectInput}
              onSubmit={handleTextAnalysis}
              isProcessing={isProcessing}
              error={error}
            />

            <ModeSelector
              cameraMode={cameraMode}
              setCameraMode={setCameraMode}
              photoMode={photoMode}
              setPhotoMode={setPhotoMode}
            />

            {cameraMode && (
              <CameraSection
                isProcessing={isProcessing}
                setIsProcessing={setIsProcessing}
                setCurrentAnalysisType={setCurrentAnalysisType}
                onAnalysisComplete={handleAnalysisComplete}
              />
            )}

            {photoMode && (
              <PhotoSection
                isProcessing={isProcessing}
                setIsProcessing={setIsProcessing}
                setCurrentAnalysisType={setCurrentAnalysisType}
                onAnalysisComplete={handleAnalysisComplete}
              />
            )}

            <Results viewers={viewers} setViewers={setViewers} />
          </div>

          {/* Right side: Payment section */}
          <PaymentSection />
        </div>
      </div>

      {/* Help button */}
      <button 
        className="help-button"
        onClick={() => setShowShortcuts(true)}
        title={`Keyboard shortcuts (${navigator.platform.toUpperCase().indexOf('MAC') >= 0 ? '⌘' : 'Ctrl'}+?)`}
      >
        ?
      </button>

      {/* Retry button for failed analyses */}
      {error && lastSuccessfulAnalysis && (
        <button 
          className="retry-button"
          onClick={handleRetry}
          title="Retry last analysis"
        >
          ↻
        </button>
      )}

      {/* Keyboard shortcuts help overlay */}
      {showShortcuts && (
        <div className="shortcuts-overlay" onClick={() => setShowShortcuts(false)}>
          <div className="shortcuts-modal" onClick={e => e.stopPropagation()}>
            <div className="shortcuts-header">
              <h3>Keyboard Shortcuts</h3>
              <button onClick={() => setShowShortcuts(false)}>×</button>
            </div>
            <div className="shortcuts-list">
              <div className="shortcut-item">
                <kbd>{navigator.platform.toUpperCase().indexOf('MAC') >= 0 ? '⌘' : 'Ctrl'}+K</kbd>
                <span>Focus text input</span>
              </div>
              <div className="shortcut-item">
                <kbd>{navigator.platform.toUpperCase().indexOf('MAC') >= 0 ? '⌘' : 'Ctrl'}+C</kbd>
                <span>Toggle camera mode</span>
              </div>
              <div className="shortcut-item">
                <kbd>{navigator.platform.toUpperCase().indexOf('MAC') >= 0 ? '⌘' : 'Ctrl'}+P</kbd>
                <span>Toggle photo mode</span>
              </div>
              <div className="shortcut-item">
                <kbd>{navigator.platform.toUpperCase().indexOf('MAC') >= 0 ? '⌘' : 'Ctrl'}+Enter</kbd>
                <span>Submit text analysis</span>
              </div>
              <div className="shortcut-item">
                <kbd>{navigator.platform.toUpperCase().indexOf('MAC') >= 0 ? '⌘' : 'Ctrl'}+⌫</kbd>
                <span>Clear all results</span>
              </div>
              <div className="shortcut-item">
                <kbd>Esc</kbd>
                <span>Clear modes & focus input</span>
              </div>
              <div className="shortcut-item">
                <kbd>{navigator.platform.toUpperCase().indexOf('MAC') >= 0 ? '⌘' : 'Ctrl'}+W</kbd>
                <span>Close last molecule viewer</span>
              </div>
              <div className="shortcut-item">
                <kbd>{navigator.platform.toUpperCase().indexOf('MAC') >= 0 ? '⌘' : 'Ctrl'}+?</kbd>
                <span>Show/hide shortcuts</span>
              </div>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default MainLayout;