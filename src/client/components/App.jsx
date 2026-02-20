import React, { useState, useEffect, useCallback } from 'react';
import { useApi } from '../hooks/useApi';
import { PRESET_VISUAL_TESTS, SMILES_NAME_MAP } from '../constants.js';
import { APP_CONFIG } from '../utils/config-loader.js';
import logger from '../logger.js';
import { createKeyboardHandler } from '../keyboard-shortcuts.js';
import ErrorBanner from './ErrorBanner';
import TextInput from './TextInput';
import CameraSection from './CameraSection';
import PhotoSection from './PhotoSection';
import LinkSection from './LinkSection';
import MoleculeViewer from './MoleculeViewer';
import MolecularColumn from './MolecularColumn';
import '../assets/style.css';

const isMobileDevice = () => window.matchMedia('(pointer: coarse) and (hover: none)').matches;
const isMac = /mac/i.test(navigator.userAgent);

const ModeSelector = ({ cameraMode, setCameraMode, photoMode, setPhotoMode, linkMode, setLinkMode }) => {
  const isMobile = isMobileDevice();

  const handleModeSelect = (mode) => {
    setCameraMode(mode === 'camera');
    setPhotoMode(mode === 'photo');
    if (setLinkMode) setLinkMode(mode === 'link');
  };

  return (
    <div className="mode-row">
      <button
        className={`mode-btn${cameraMode ? ' active' : ''}`}
        onClick={() => handleModeSelect('camera')}
        title={isMobile ? "Capture from camera" : "Capture from camera (⌘⇧C)"}
      >
        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
          <circle cx="12" cy="12" r="5" fill="currentColor" opacity="0.8"/>
          <circle cx="12" cy="12" r="9" stroke="currentColor" fill="none"/>
        </svg>
        {isMobile && <span className="mode-label">Live</span>}
        {!isMobile && <span className="mode-btn-shortcut">{isMac ? '⌘⇧C' : 'Ctrl+Shift+C'}</span>}
      </button>
      <button
        className={`mode-btn${photoMode ? ' active' : ''}`}
        onClick={() => handleModeSelect('photo')}
        title={isMobile ? "Upload image" : "Upload image (⌘⇧P)"}
      >
        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
          <rect x="3" y="3" width="18" height="18" rx="2" ry="2"/>
          <circle cx="9" cy="9" r="2"/>
          <path d="M21 15l-3.086-3.086a2 2 0 0 0-2.828 0L6 21"/>
        </svg>
        {isMobile && <span className="mode-label">Gallery</span>}
        {!isMobile && <span className="mode-btn-shortcut">{isMac ? '⌘⇧P' : 'Ctrl+Shift+P'}</span>}
      </button>
      <button
        className={`mode-btn${linkMode ? ' active' : ''}`}
        onClick={() => handleModeSelect('link')}
        title={isMobile ? "Enter image link" : "Enter image link (⌥L)"}
      >
        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
          <path d="M10 13a5 5 0 0 0 7.54.54l3-3a5 5 0 0 0-7.07-7.07l-1.72 1.71"/>
          <path d="M14 11a5 5 0 0 0-7.54-.54l-3 3a5 5 0 0 0 7.07 7.07l1.71-1.71"/>
        </svg>
        {isMobile && <span className="mode-label">Link</span>}
        {!isMobile && <span className="mode-btn-shortcut">{isMac ? '⌘⇧L' : 'Ctrl+Shift+L'}</span>}
      </button>
    </div>
  );
};

function App() {
  const [objectInput, setObjectInput] = useState(APP_CONFIG.defaultInput || '');
  const [cameraMode, setCameraMode] = useState(false);
  const [photoMode, setPhotoMode] = useState(false);
  const [linkMode, setLinkMode] = useState(false);
  const [isProcessing, setIsProcessing] = useState(false);
  const [columns, setColumns] = useState([]);
  const [error, setError] = useState('');
  // Visual tests only run when explicitly enabled (not automatic in dev mode)
  // Enable via: REACT_APP_RUN_VISUAL_TESTS=true or URL ?visual-tests
  const [autoVisualMode, setAutoVisualMode] = useState(
    process.env.NODE_ENV === 'development' &&
    (process.env.REACT_APP_RUN_VISUAL_TESTS === 'true' || new URLSearchParams(window.location.search).has('visual-tests'))
  );

  // Check for test data loading via URL parameters
  const urlParams = new URLSearchParams(window.location.search);
  const testName = urlParams.get('test');
  const autoload = urlParams.get('autoload') === 'true';
  const [showSettings, setShowSettings] = useState(false);
  const [columnMode, setColumnMode] = useState('accumulate'); // 'replace' or 'accumulate'
  const [lookupMode, setLookupMode] = useState('GPT-5');
  const [showLeftSidebar, setShowLeftSidebar] = useState(true); // New state for sidebar visibility

  const { structuresFromText: analyzeText, generateSDFs } = useApi();

  // Load 3Dmol.js once at app start
  useEffect(() => {
    if (typeof window.$3Dmol !== 'undefined') return;
    if (document.getElementById('threedmol-script')) return;
    const script = document.createElement('script');
    script.id = 'threedmol-script';
    script.src = 'https://3Dmol.org/build/3Dmol-min.js';
    script.async = true;
    document.head.appendChild(script);
  }, []);

  // Auto-enable dev mode
  useEffect(() => {
    // Dev-mode setup only
  }, []);

  // Ensure accumulate mode during visual tests
  useEffect(() => {
    if (autoVisualMode && columnMode !== 'accumulate') {
      setColumnMode('accumulate');
    }
  }, [autoVisualMode, columnMode]);

  // Load test data if specified in URL parameters
  useEffect(() => {
    if (testName && autoload) {
      const loadTestData = async () => {
        try {
          setIsProcessing(true);
          setError('');

          // Fetch test data from server
          const response = await fetch(`/api/test-data/${testName}`);
          if (!response.ok) {
            throw new Error(`HTTP ${response.status}: ${response.statusText} (endpoint: /api/test-data/${testName})`);
          }

          const testData = await response.json();

          // Set the object input to show what was tested
          setObjectInput(testData.object || '');

          // Create visualization data from test results
          const visualizationData = testData.chemicals?.map(chem => ({
            name: chem.name,
            smiles: null, // We'll use SDF files instead
            sdfData: chem.sdfUrl ? `file://${chem.sdfUrl}` : null,
            formula: null, // Could be enhanced later
            status: chem.status
          })).filter(mol => mol.sdfData) || [];

          // Create a column with the test results
          const newColumn = {
            id: `test-${testName}-${Date.now()}`,
            title: `Test: ${testData.object}`,
            subtitle: `Generated: ${new Date(testData.timestamp).toLocaleString()}`,
            objectName: testData.object,
            molecules: visualizationData,
            metadata: testData.metadata,
            isTestData: true
          };

          // Set columns based on column mode
          if (columnMode === 'replace') {
            setColumns([newColumn]);
          } else {
            setColumns(prev => [...prev, newColumn]);
          }

        } catch (error) {
          const isDev = process.env.NODE_ENV === 'development';
          const errorMsg = `Test data load failed: ${error.message || error.toString()} (test: ${testName}, status: ${error.response?.status || 'N/A'})`;
          if (isDev) console.error('[App] Test data error', { error, testName, response: error.response });
          setError(errorMsg);
        } finally {
          setIsProcessing(false);
        }
      };

      loadTestData();
    }
  }, [testName, autoload, columnMode]); // Re-run if URL parameters change

  // Handle keyboard shortcuts (desktop only)
  useEffect(() => {
    // Skip keyboard shortcuts on mobile devices
    if (isMobileDevice()) {
      logger.info('Skipping shortcuts - mobile device detected');
      return;
    }

    // Create keyboard handler with mode actions
    const actions = {
      focusInput: () => {
        const inputElement = document.getElementById('object-input');
        if (inputElement) {
          inputElement.focus();
          inputElement.select();
        }
      },
      cameraMode: () => {
        setCameraMode(true);
        setPhotoMode(false);
        setLinkMode(false);
      },
      photoMode: () => {
        setCameraMode(false);
        setPhotoMode(true);
        setLinkMode(false);
      },
      linkMode: () => {
        setCameraMode(false);
        setPhotoMode(false);
        setLinkMode(true);
      }
    };

    const keyboardHandler = createKeyboardHandler(actions);

    document.addEventListener('keydown', keyboardHandler, { capture: true });

    return () => document.removeEventListener('keydown', keyboardHandler, { capture: true });
  }, [setCameraMode, setPhotoMode, setLinkMode]);

  const updateColumn = useCallback((columnId, updates) => {
    setColumns(prev => prev.map(col => col.id === columnId ? { ...col, ...updates } : col));
  }, []);

  const handleTextPrediction = useCallback(async (value) => {
    setIsProcessing(true);
    setError('');

    const columnId = Date.now();
    const newColumn = { id: columnId, query: value, viewers: [], loading: true, failed: false };

    setColumns(prev => columnMode === 'replace' ? [newColumn] : [...prev, newColumn]);

    try {
      const result = await analyzeText(value, lookupMode);
      const molecules = result.molecules || result.chemicals || [];

      if (molecules && molecules.length > 0) {
        // Prefer precomputed SDFs with canonical names
        const precomputed = molecules.filter(m => m.sdfPath);
        if (precomputed.length > 0) {
          const viewers = precomputed.map(m => ({
            name: m.name || value,
            sdfData: m.sdfPath.startsWith('file://') ? m.sdfPath : `file://${m.sdfPath}`,
            smiles: m.smiles
          }));
          updateColumn(columnId, { viewers, loading: false, failed: false });
        } else {
          // Legacy path: generate SDFs from SMILES
          const smilesArray = molecules.map(mol => mol.smiles).filter(Boolean);
          if (smilesArray.length > 0) {
            try {
              const sdfResult = await generateSDFs(smilesArray, false);
              const smilesToSdf = new Map();
              sdfResult.sdfPaths?.forEach((p, i) => { smilesToSdf.set(smilesArray[i], p); });
              const viewers = molecules.map((mol) => {
                const sdfPath = smilesToSdf.get(mol.smiles);
                return {
                  name: mol.name || (mol.smiles && SMILES_NAME_MAP[mol.smiles]) || mol.smiles || value,
                  sdfData: sdfPath ? `file://${sdfPath}` : null,
                  smiles: mol.smiles
                };
              });
              updateColumn(columnId, { viewers, loading: false, failed: false });
            } catch (sdfError) {
              const isDev = process.env.NODE_ENV === 'development';
              const errorDetails = sdfError.details || {};
              const errorMsg = `generateSDFs() timed out: ${sdfError.message || sdfError.toString() || 'Unknown error'} (endpoint: /api/generate-sdfs, SMILES count: ${smilesArray.length}, input: "${value}")`;
              if (isDev) console.error('[handleTextPrediction] SDF generation error', {
                sdfError,
                errorDetails,
                smilesArray,
                columnId,
                function: 'generateSDFs',
                endpoint: '/api/generate-sdfs'
              });
              logger.error('SDF generation failed:', sdfError);
              setError(errorMsg);
              updateColumn(columnId, { loading: false, failed: true });
            }
          } else {
            const isDev = process.env.NODE_ENV === 'development';
            const errorMsg = `No SMILES found in ${molecules.length} molecules from API`;
            setError(errorMsg);
            updateColumn(columnId, { loading: false, failed: true });
          }
        }
      } else {
        updateColumn(columnId, { viewers: [], loading: false, failed: false });
      }

      setObjectInput('');
      } catch (error) {
      const isDev = process.env.NODE_ENV === 'development';
      const errorDetails = error.details || {};
      const isTimeout = error.message?.includes('timeout') || error.message?.includes('timed out') || errorDetails.message?.includes('timeout');
      const failedFunction = errorDetails.endpoint?.includes('generate-sdfs') ? 'generateSDFs()' : 'analyzeText()';
      const failedEndpoint = errorDetails.endpoint || '/api/structuralize';

      if (isDev) {
        console.error('[handleTextPrediction] ERROR', {
          error,
          errorDetails,
          value,
          lookupMode,
          columnId,
          failedFunction,
          failedEndpoint,
          stack: error.stack
        });
      }
      logger.error('Prediction failed:', error);

      let errorMessage;
      if (isTimeout) {
        errorMessage = `${failedFunction} timed out: ${error.message || errorDetails.message || 'Request timed out'} (endpoint: ${failedEndpoint}, input: "${value}", mode: ${lookupMode})`;
      } else if (error.details?.message) {
        errorMessage = `${failedFunction} failed: ${error.details.message} (endpoint: ${failedEndpoint}, input: "${value}")`;
      } else if (error.message) {
        errorMessage = `${failedFunction} failed: ${error.message} (endpoint: ${failedEndpoint}, input: "${value}")`;
      } else if (error.toString) {
        errorMessage = `${failedFunction} error: ${error.toString()} (endpoint: ${failedEndpoint}, input: "${value}")`;
      } else {
        errorMessage = `${failedFunction} unknown error: ${JSON.stringify(error)} (endpoint: ${failedEndpoint}, input: "${value}")`;
      }

      setError(errorMessage);
      updateColumn(columnId, { loading: false, failed: true });
    } finally {
      setIsProcessing(false);
    }
  }, [analyzeText, generateSDFs, columnMode, updateColumn]);

  const handlePredictionComplete = useCallback(async (result) => {
    const molecules = result?.molecules || result?.chemicals || [];
    const objectLabel = result?.object || (cameraMode ? 'Camera capture' : 'Image capture');
    const columnId = Date.now();
    const newColumn = { id: columnId, query: objectLabel, viewers: [], loading: true, failed: false };

    setColumns(prev => columnMode === 'replace' ? [newColumn] : [...prev, newColumn]);

    if (molecules?.length > 0) {
      const smilesArray = molecules.map(mol => mol.smiles).filter(Boolean);
      if (smilesArray.length > 0) {
        try {
          const sdfResult = await generateSDFs(smilesArray, false);
          const smilesToSdf = new Map();
          sdfResult.sdfPaths?.forEach((p, i) => { smilesToSdf.set(smilesArray[i], p); });
          const viewers = molecules.map((mol) => ({
            name: mol.name || (mol.smiles && SMILES_NAME_MAP[mol.smiles]) || mol.smiles || objectLabel,
            sdfData: smilesToSdf.get(mol.smiles) ? `file://${smilesToSdf.get(mol.smiles)}` : null,
            smiles: mol.smiles
          }));
          updateColumn(columnId, { viewers, loading: false, failed: false });
        } catch (sdfError) {
          const isDev = process.env.NODE_ENV === 'development';
          const errorMsg = `SDF generation failed: ${sdfError.message || sdfError.toString() || 'Unknown error'}`;
          if (isDev) console.error('[handlePredictionComplete] SDF error', { sdfError, smilesArray, columnId });
          logger.error('SDF generation failed:', sdfError);
          setError(errorMsg);
          updateColumn(columnId, { loading: false, failed: true });
        }
      } else {
        updateColumn(columnId, { loading: false, failed: true });
      }
    } else {
      updateColumn(columnId, { loading: false, failed: true });
    }
    setError('');
  }, [generateSDFs, cameraMode, columnMode, updateColumn]);

  // Auto-run visual tests - bypass AI and go straight to SDF generation
  useEffect(() => {
    let cancelled = false;
    const run = async () => {
      if (!autoVisualMode) return;

      // Wait a bit to ensure 3Dmol.js is loaded and DOM is ready
      await new Promise(resolve => setTimeout(resolve, 1000));

      for (const test of PRESET_VISUAL_TESTS) {
        if (cancelled) break;
        try {
          // Generate SDFs for multiple compounds per object
          const smilesArray = test.smilesList || [];
          if (smilesArray.length === 0) continue;
          const sdfResult = await generateSDFs(smilesArray, false);
          // Prefer actual returned paths; fall back to deterministic filenames
          const sanitizeSmiles = (s) => s.replace(/[^a-zA-Z0-9]/g, ch => ch === '=' ? '__' : '_');
          const returnedPaths = Array.isArray(sdfResult?.sdfPaths) ? sdfResult.sdfPaths : [];
          const returnedFileSet = new Set(returnedPaths.map(p => p.replace(/^\/*/, '/')));
          const viewers = smilesArray.map((sm, idx) => {
            const deterministic = `/sdf_files/${sanitizeSmiles(sm)}.sdf`;
            const byIndex = returnedPaths[idx] || null;
            const chosen = returnedFileSet.has(deterministic) ? deterministic : (byIndex || deterministic);
            return {
              name: SMILES_NAME_MAP[sm] || sm,
              sdfData: chosen ? `file://${chosen}` : null,
              smiles: sm
            };
          });
          setColumns(prev => ([...prev, { id: Date.now() + Math.random(), query: test.label, viewers }]));
        } catch (e) {}
      }
    };
    run();
    return () => { cancelled = true; };
  }, [autoVisualMode, generateSDFs]);

  const removeColumn = useCallback((columnId) => {
    setColumns(prev => prev.filter(col => col.id !== columnId));
  }, []);



  return (
    <div className="app">
      {/* Desktop split-screen layout */}
      {!isMobileDevice() ? (
        <div className="split-screen-container">
          {/* Toggle button for left sidebar */}
          <button
            className="sidebar-toggle"
            onClick={() => setShowLeftSidebar(!showLeftSidebar)}
            title={showLeftSidebar ? 'Hide sidebar' : 'Show sidebar'}
            style={{ left: showLeftSidebar ? '400px' : '0' }}
          >
            {showLeftSidebar ? '◀' : '▶'}
          </button>

          {/* Left sidebar with input modes */}
          <div className={`left-sidebar ${showLeftSidebar ? 'visible' : 'hidden'}`}>

            <div className="input-section">
                <TextInput
                  value={objectInput}
                  onChange={setObjectInput}
                  onSubmit={handleTextPrediction}
                  isProcessing={isProcessing}
                  error={error}
                />

                <ModeSelector
                  cameraMode={cameraMode}
                  setCameraMode={setCameraMode}
                  photoMode={photoMode}
                  setPhotoMode={setPhotoMode}
                  linkMode={linkMode}
                  setLinkMode={setLinkMode}
                  lookupMode={lookupMode}
                  setLookupMode={setLookupMode}
                />

                {cameraMode && (
                  <CameraSection
                    isProcessing={isProcessing}
                    setIsProcessing={setIsProcessing}
                    setCurrentAnalysisType={() => {}}
                    onAnalysisComplete={handlePredictionComplete}
                  />
                )}

                {photoMode && (
                  <PhotoSection
                    isProcessing={isProcessing}
                    setIsProcessing={setIsProcessing}
                    setCurrentAnalysisType={() => {}}
                    onAnalysisComplete={handlePredictionComplete}
                  />
                )}

                {linkMode && (
                  <LinkSection
                    isProcessing={isProcessing}
                    setIsProcessing={setIsProcessing}
                    onPredictionComplete={handlePredictionComplete}
                  />
                )}
              </div>
          </div>

          {/* Right side with molecule viewer and settings */}
          <div className="right-content">
            <ErrorBanner error={error} onDismiss={() => setError('')} />

            {/* Settings gear icon in top right */}
            <button
              onClick={() => {
                logger.debug('Settings button clicked', { showSettings });
                setShowSettings(!showSettings);
              }}
              className="settings-btn"
              title="Settings"
            >
              ⚙️
            </button>

            {/* Settings Modal */}
            {showSettings && (
              <div className="settings-modal">
                <h3 className="settings-title">Settings</h3>

                <div className="settings-field">
                  <label className="label">
                    Column Behavior:
                  </label>
                  <select
                    value={columnMode}
                    onChange={(e) => {
                      logger.info('Column mode changed', { value: e.target.value });
                      setColumnMode(e.target.value);
                    }}
                    className="select"
                  >
                    <option value="replace">
                      Replace (Default) - New prediction displaces previous
                    </option>
                    <option value="accumulate">
                      Accumulate - Add columns side by side
                    </option>
                  </select>
                </div>

                {process.env.NODE_ENV === 'development' && (
                  <div className="settings-field">
                    <label className="label">
                      <input
                        type="checkbox"
                        checked={autoVisualMode}
                        onChange={(e) => setAutoVisualMode(e.target.checked)}
                        className="checkbox-input"
                      />
                      Enable Visual Tests (Dev Mode)
                    </label>
                  </div>
                )}

                <button
                  onClick={() => setShowSettings(false)}
                  className="settings-close"
                >
                  Close
                </button>
              </div>
            )}

            {/* Molecule viewer columns */}
            <div className="columns">
              {columns.map(column => (
                <MolecularColumn
                  key={column.id}
                  column={column}
                  onRemove={() => removeColumn(column.id)}
                  showRemove={columnMode === 'accumulate'}
                />
              ))}
            </div>
          </div>
        </div>
      ) : (
        /* Mobile layout remains unchanged */
        <div className="main mobile-reordered">
          <h1 className="app-title">what's in...</h1>

          <ErrorBanner error={error} onDismiss={() => setError('')} />

          <div className="columns mobile-top">
            {columns.map(column => (
              <MolecularColumn
                key={column.id}
                column={column}
                onRemove={() => removeColumn(column.id)}
                showRemove={columnMode === 'accumulate'}
              />
            ))}
          </div>

          <div className="input-section mobile-bottom">
            <TextInput
              value={objectInput}
              onChange={setObjectInput}
              onSubmit={handleTextPrediction}
              isProcessing={isProcessing}
              error={error}
            />

            <ModeSelector
              cameraMode={cameraMode}
              setCameraMode={setCameraMode}
              photoMode={photoMode}
              setPhotoMode={setPhotoMode}
              linkMode={linkMode}
              setLinkMode={setLinkMode}
              lookupMode={lookupMode}
              setLookupMode={setLookupMode}
            />

            {cameraMode && (
              <CameraSection
                isProcessing={isProcessing}
                setIsProcessing={setIsProcessing}
                setCurrentAnalysisType={() => {}}
                onAnalysisComplete={handlePredictionComplete}
              />
            )}

            {photoMode && (
              <PhotoSection
                isProcessing={isProcessing}
                setIsProcessing={setIsProcessing}
                setCurrentAnalysisType={() => {}}
                onAnalysisComplete={handlePredictionComplete}
              />
            )}

            {linkMode && (
              <LinkSection
                isProcessing={isProcessing}
                setIsProcessing={setIsProcessing}
                onPredictionComplete={handlePredictionComplete}
              />
            )}
          </div>
        </div>
      )}
    </div>
  );
}

export default App;
