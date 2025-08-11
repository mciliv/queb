import React, { useState, useEffect, useCallback, useRef } from 'react';
import TextInput from '../input/TextInput';
import ModeSelector from '../input/ModeSelector';
import CameraSection from '../input/CameraSection';
import PhotoSection from '../input/PhotoSection';
import LinkSection from '../input/LinkSection';
import { useApi } from '../../hooks/useApi';

// Inline styles for main layout
const styles = {
  appContainer: {
    background: '#000000',
    color: '#ffffff',
    minHeight: '100vh',
    width: '100vw',
    maxWidth: '100vw',
    overflowX: 'clip',
    fontFamily: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif',
    fontSize: '13px',
    fontWeight: 400,
    lineHeight: 1.4,
    letterSpacing: '0.01em'
  },
  mainLayout: {
    display: 'flex',
    flexDirection: 'column',
    padding: '20px 0',
    gap: '20px'
  },
  visualToggleButton: {
    position: 'fixed',
    left: '12px',
    bottom: '12px',
    width: '40px',
    height: '40px',
    borderRadius: '20px',
    background: 'rgba(255, 255, 255, 0.08)',
    color: '#ffffff',
    border: 'none',
    outline: 'none',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    cursor: 'pointer'
  },
  inputSection: {
    display: 'flex',
    flexDirection: 'column',
    gap: '10px',
    maxWidth: '600px',
    paddingLeft: '20px',
    paddingRight: '20px'
  },
  columnsContainer: {
    display: 'flex',
    flexDirection: 'row',
    gap: '20px',
    overflowX: 'auto',
    padding: '20px',
    minHeight: '400px'
  },
  column: {
    minWidth: '400px',
    background: 'transparent',
    padding: '20px'
  }
};

// Preset visual tests (module-scoped constant, no refs)
const PRESET_VISUAL_TESTS = [
  { label: 'Water', smiles: 'O' },
  { label: 'Ethanol', smiles: 'CCO' },
  { label: 'Acetic acid', smiles: 'CC(=O)O' }
];

const MainLayout = () => {
  const [objectInput, setObjectInput] = useState('');
  const [cameraMode, setCameraMode] = useState(false);
  const [photoMode, setPhotoMode] = useState(false);
  const [linkMode, setLinkMode] = useState(false);
  const [isProcessing, setIsProcessing] = useState(false);
  const [columns, setColumns] = useState([]);
  const [error, setError] = useState('');
  const [autoVisualMode, setAutoVisualMode] = useState(true);
  const { analyzeText, generateSDFs } = useApi();

  // Load 3Dmol.js once at app start (simplified)
  useEffect(() => {
    if (typeof window.$3Dmol !== 'undefined') return;
    if (document.getElementById('threedmol-script')) return;
    const script = document.createElement('script');
    script.id = 'threedmol-script';
    script.src = 'https://3Dmol.org/build/3Dmol-min.js';
    script.async = true;
    document.head.appendChild(script);
  }, []);

  // Handle keyboard shortcuts
  useEffect(() => {
    const handleKeyDown = (event) => {
      if (event.target.tagName === 'INPUT' || event.target.tagName === 'TEXTAREA') {
        return;
      }

      const isMac = navigator.platform.toUpperCase().indexOf('MAC') >= 0;
      const modifier = isMac ? event.metaKey : event.ctrlKey;

      if (modifier && event.key.toLowerCase() === 'k') {
        event.preventDefault();
        document.getElementById('object-input')?.focus();
      }
    };

    document.addEventListener('keydown', handleKeyDown);
    return () => document.removeEventListener('keydown', handleKeyDown);
  }, []);

  const handleTextAnalysis = useCallback(async (value) => {
    if (!value.trim()) return;

    setIsProcessing(true);
    setError('');
    
    try {
      const result = await analyzeText(value);
      const molecules = result.molecules || result.chemicals || [];
      
      if (molecules && molecules.length > 0) {
        const smilesArray = molecules.map(mol => mol.smiles).filter(Boolean);
        
        if (smilesArray.length > 0) {
          try {
            const sdfResult = await generateSDFs(smilesArray, false);
            
            const viewers = molecules.map((mol, index) => {
              const sdfPath = sdfResult.sdfPaths && sdfResult.sdfPaths[index];
              return {
                name: mol.name || value,
                sdfData: sdfPath ? `file://${sdfPath}` : null,
                smiles: mol.smiles
              };
            });
            
            // Add molecules to existing column or create new one
            setColumns(prev => {
              if (prev.length === 0) {
                // First analysis - create new column
                return [{
                  id: Date.now(),
                  query: value,
                  viewers: viewers
                }];
              } else {
                // Subsequent analyses - add to existing column
                const updatedColumns = [...prev];
                const lastColumn = updatedColumns[updatedColumns.length - 1];
                lastColumn.viewers = [...lastColumn.viewers, ...viewers];
                lastColumn.query = `${lastColumn.query} + ${value}`;
                return updatedColumns;
              }
            });
          } catch (sdfError) {
            console.error('SDF generation failed:', sdfError);
            setError('Failed to generate molecular structures');
          }
        } else {
          setError('No valid SMILES found in the analysis results.');
        }
      } else {
        setError('No molecules found for this input.');
      }
      
      setObjectInput('');
    } catch (error) {
      console.error('Analysis failed:', error);
      setError('Analysis failed. Please try again.');
    } finally {
      setIsProcessing(false);
    }
  }, [isProcessing, analyzeText, generateSDFs]);

  const handleAnalysisComplete = useCallback(async (result) => {
    const molecules = result.molecules || result.chemicals || [];
    
    if (molecules && molecules.length > 0) {
      const smilesArray = molecules.map(mol => mol.smiles).filter(Boolean);
      
      if (smilesArray.length > 0) {
        try {
          const sdfResult = await generateSDFs(smilesArray, false);
          
          const viewers = molecules.map((mol, index) => {
            const sdfPath = sdfResult.sdfPaths && sdfResult.sdfPaths[index];
            return {
              name: mol.name || 'Captured object',
              sdfData: sdfPath ? `file://${sdfPath}` : null,
              smiles: mol.smiles
            };
          });
          
          // Add molecules to existing column or create new one  
          setColumns(prev => {
            const captureType = cameraMode ? 'camera' : 'image';
            if (prev.length === 0) {
              // First analysis - create new column
              return [{
                id: Date.now(),
                query: `Captured from ${captureType}`,
                viewers: viewers
              }];
            } else {
              // Subsequent analyses - add to existing column
              const updatedColumns = [...prev];
              const lastColumn = updatedColumns[updatedColumns.length - 1];
              lastColumn.viewers = [...lastColumn.viewers, ...viewers];
              lastColumn.query = `${lastColumn.query} + Captured from ${captureType}`;
              return updatedColumns;
            }
          });
        } catch (sdfError) {
          console.error('SDF generation failed:', sdfError);
        }
      }
    }
    
    setError('');
  }, [generateSDFs, cameraMode]);

  // Predefined visual test inputs (module constant)

  // Auto-run visual tests sequentially, appending each as a column
  useEffect(() => {
    let cancelled = false;
    const run = async () => {
      if (!autoVisualMode) return;
      // If any columns already exist from manual usage, keep appending a new column per test
      for (const test of PRESET_VISUAL_TESTS) {
        if (cancelled) break;
        try {
          // Analyze text path to reuse existing pipeline and SDF generation
          const result = await analyzeText(test.smiles);
          const molecules = result.molecules || result.chemicals || [];
          const smilesArray = molecules.map(m => m.smiles).filter(Boolean);
          if (smilesArray.length > 0) {
            const sdfResult = await generateSDFs(smilesArray, false);
            const viewers = molecules.map((mol, index) => {
              const sdfPath = sdfResult.sdfPaths && sdfResult.sdfPaths[index];
              return {
                name: mol.name || test.label,
                sdfData: sdfPath ? `file://${sdfPath}` : null,
                smiles: mol.smiles
              };
            });
            setColumns(prev => ([...prev, { id: Date.now() + Math.random(), query: test.label, viewers }]));
          }
        } catch (e) {
          // Non-fatal; continue
        }
      }
    };
    run();
    return () => { cancelled = true; };
  }, [autoVisualMode, analyzeText, generateSDFs]);

  const removeColumn = useCallback((columnId) => {
    setColumns(prev => prev.filter(col => col.id !== columnId));
  }, []);

  return (
    <div style={styles.appContainer}>
      <div style={styles.mainLayout}>
        {/* Visual tests toggle (beaker icon) */}
        <button
          aria-label="Toggle visual tests"
          title="Toggle visual tests"
          onClick={() => setAutoVisualMode(v => !v)}
          style={styles.visualToggleButton}
        >
          {/* Beaker unicode */}
          ⚗️
        </button>
        {/* Input section */}
        <div style={styles.inputSection}>
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
            linkMode={linkMode}
            setLinkMode={setLinkMode}
          />

          {cameraMode && (
            <CameraSection
              isProcessing={isProcessing}
              setIsProcessing={setIsProcessing}
              setCurrentAnalysisType={() => {}}
              onAnalysisComplete={handleAnalysisComplete}
            />
          )}

          {photoMode && (
            <PhotoSection
              isProcessing={isProcessing}
              setIsProcessing={setIsProcessing}
              setCurrentAnalysisType={() => {}}
              onAnalysisComplete={handleAnalysisComplete}
            />
          )}

          {linkMode && (
            <LinkSection
              isProcessing={isProcessing}
              setIsProcessing={setIsProcessing}
              onAnalysisComplete={handleAnalysisComplete}
            />
          )}
        </div>

        {/* Columns container */}
        <div style={styles.columnsContainer}>
          {columns.map(column => (
            <MolecularColumn
              key={column.id}
              column={column}
              onRemove={() => removeColumn(column.id)}
            />
          ))}
        </div>
      </div>
    </div>
  );
};

// Minimal per-molecule viewer component
const SimpleMoleculeViewer = ({ molecularData }) => {
  const ref = useRef(null);

  useEffect(() => {
    let cancelled = false;
    const initialize = async () => {
      while (typeof window.$3Dmol === 'undefined') {
        await new Promise(resolve => setTimeout(resolve, 100));
        if (cancelled) return;
      }
      if (!ref.current) return;

      const viewer = window.$3Dmol.createViewer(ref.current, {
        backgroundColor: 'transparent',
        antialias: true,
        defaultcolors: window.$3Dmol.rasmolElementColors
      });

      let sdfContent = null;
      if (molecularData.sdfData && molecularData.sdfData.startsWith('file://')) {
        const path = molecularData.sdfData.replace('file://', '');
        const response = await fetch(`/sdf_files/${path.split('/').pop()}`);
        if (response.ok) {
          sdfContent = await response.text();
        }
      }

      if (sdfContent) {
        viewer.addModel(sdfContent, 'sdf');
        viewer.setStyle({}, { sphere: { scale: 0.8 } });
        viewer.zoomTo();
        viewer.render();
      }
    };

    initialize();
    return () => {
      cancelled = true;
      if (ref.current) {
        ref.current.innerHTML = '';
      }
    };
  }, [molecularData.sdfData]);

  return (
    <div style={{ marginBottom: '12px' }}>
      <div 
        ref={ref}
        style={{ height: '200px', width: '100%', background: 'transparent', border: 'none' }}
      />
    </div>
  );
};

// Individual column component (simplified)
const MolecularColumn = ({ column, onRemove }) => {
  return (
    <div style={styles.column}>
      <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '20px' }}>
        <div aria-hidden role="img" style={{ opacity: 0.8 }}>🧪</div>
        <button 
          onClick={onRemove}
          style={{
            background: 'transparent',
            border: 'none',
            color: '#ffffff',
            fontSize: '20px',
            cursor: 'pointer'
          }}
        >
          ×
        </button>
      </div>

      {column.viewers.map((mol, idx) => (
        <SimpleMoleculeViewer key={idx} molecularData={mol} />
      ))}
    </div>
  );
};

export default MainLayout;