import React, { useState, useEffect, useRef, useCallback } from 'react';
import TextInput from '../input/TextInput';
import ModeSelector from '../input/ModeSelector';
import CameraSection from '../input/CameraSection';
import PhotoSection from '../input/PhotoSection';
import MolecularAnalysisResults from '../visualization/Results';
import PaymentSection from './PaymentSection';
import MolecularTestPanel from '../visualization/MolecularTestPanel';
import { usePayment } from './PaymentContext';
import { useApi } from '../../hooks/useApi';

// Inline styles for main layout
const styles = {
  appContainer: {
    background: '#000000',
    color: '#ffffff',
    minHeight: '100vh',
    fontFamily: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif',
    fontSize: '13px',
    fontWeight: 400,
    lineHeight: 1.4,
    letterSpacing: '0.01em'
  },
  mainAppInterface: {
    position: 'relative',
    width: '100%',
    minHeight: '100vh'
  },
  mainContentLayout: {
    display: 'flex',
    flexDirection: 'row',
    gap: '20px',
    padding: '20px',
    maxWidth: '100vw',
    overflow: 'hidden'
  },
  analysisSection: {
    flex: 1,
    display: 'flex',
    flexDirection: 'column',
    gap: '20px',
    maxWidth: 'calc(100vw - 320px)'
  },
  helpButton: {
    position: 'fixed',
    bottom: '20px',
    left: '20px',
    background: 'rgba(255, 255, 255, 0.1)',
    border: 'none',
    borderRadius: '50%',
    width: '50px',
    height: '50px',
    color: '#ffffff',
    cursor: 'pointer',
    fontSize: '18px',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    zIndex: 1000
  }
};

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
  // API hook for connection testing
  const { testConnection } = useApi();
  
  // Debug: Test API connection on component mount
  useEffect(() => {
    const runConnectionTest = async () => {
      try {
        console.log('🔍 Testing API connection...');
        const isConnected = await testConnection();
        if (isConnected) {
          console.log('✅ API connection successful');
        } else {
          console.log('❌ API connection failed');
        }
      } catch (error) {
        console.log('❌ API connection test failed:', error.message);
      }
    };
    
    runConnectionTest();
  }, [testConnection]);

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
        case 's':
          if (modifier) {
            event.preventDefault();
            setShowShortcuts(prev => !prev);
          }
          break;
        case 't':
          if (modifier) {
            event.preventDefault();
            // Trigger test panel - find and click the test button
            const testButton = document.querySelector('.test-panel-toggle');
            if (testButton) {
              testButton.click();
            }
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
      console.log('🧪 Full API response:', result); // Debug the response
      
      // Handle both 'molecules' and 'chemicals' field names for compatibility
      const molecules = result.molecules || result.chemicals || [];
      console.log('🧪 Extracted molecules:', molecules); // Debug molecule extraction
      
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
            
            setViewers(newViewers);
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
            setViewers(newViewers);
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

  const handleRetry = useCallback(() => {
    if (lastSuccessfulAnalysis) {
      handleTextAnalysis(objectInput || lastSuccessfulAnalysis.query);
    }
  }, [lastSuccessfulAnalysis, objectInput, handleTextAnalysis]);

  const handleTestAnalysis = useCallback(async (testInput) => {
    console.log(`🧪 Running test analysis for: ${testInput}`);
    await handleTextAnalysis(testInput);
  }, [handleTextAnalysis]);

  return (
    <div style={styles.appContainer}>
      <div style={styles.mainAppInterface}>
        <div style={styles.mainContentLayout}>
          {/* Left side: Analysis section */}
          <div style={styles.analysisSection}>
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

            <MolecularAnalysisResults 
              moleculeViewers={viewers} 
              setMoleculeViewers={setViewers} 
              lastAnalysisResult={lastAnalysis}
              isAnalyzing={isProcessing}
              currentAnalysisType={currentAnalysisType}
              targetObjectInput={objectInput}
            />
          </div>

          {/* Right side: Payment section */}
          <PaymentSection />
        </div>
      </div>

      {/* Help button */}
      <button 
        style={styles.helpButton}
        onClick={() => setShowShortcuts(true)}
        title={`Keyboard shortcuts (${navigator.platform.toUpperCase().indexOf('MAC') >= 0 ? '⌘' : 'Ctrl'}+S)`}
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
                <kbd>{navigator.platform.toUpperCase().indexOf('MAC') >= 0 ? '⌘' : 'Ctrl'}+S</kbd>
                <span>Show/hide shortcuts</span>
              </div>
              <div className="shortcut-item">
                <kbd>{navigator.platform.toUpperCase().indexOf('MAC') >= 0 ? '⌘' : 'Ctrl'}+T</kbd>
                <span>Open test panel</span>
              </div>
            </div>
          </div>
        </div>
      )}

      {/* Molecular Test Panel for development */}
      <MolecularTestPanel 
        onTestAnalysis={handleTestAnalysis}
        isProcessing={isProcessing}
      />
    </div>
  );
};

export default MainLayout;