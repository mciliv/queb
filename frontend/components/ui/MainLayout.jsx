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

// Global flag to ensure 3Dmol is only loaded once
let threeDmolLoading = false;
let threeDmolLoaded = false;

const MainLayout = () => {
  const [objectInput, setObjectInput] = useState('');
  const [cameraMode, setCameraMode] = useState(false);
  const [photoMode, setPhotoMode] = useState(false);
  const [linkMode, setLinkMode] = useState(false);
  const [isProcessing, setIsProcessing] = useState(false);
  const [columns, setColumns] = useState([]);
  const [error, setError] = useState('');
  const { analyzeText, generateSDFs } = useApi();

  // Load 3Dmol.js once at app start
  useEffect(() => {
    if (!threeDmolLoaded && !threeDmolLoading && typeof window.$3Dmol === 'undefined') {
      threeDmolLoading = true;
      const script = document.createElement('script');
      script.src = 'https://3Dmol.org/build/3Dmol-min.js';
      script.async = true;
      script.onload = () => {
        threeDmolLoaded = true;
        threeDmolLoading = false;
        console.log('3Dmol.js loaded successfully');
      };
      document.head.appendChild(script);
    }
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

  const removeColumn = useCallback((columnId) => {
    setColumns(prev => prev.filter(col => col.id !== columnId));
  }, []);

  return (
    <div style={styles.appContainer}>
      <div style={styles.mainLayout}>
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

// Individual column component
const MolecularColumn = ({ column, onRemove }) => {
  const gridViewerRef = useRef(null);
  const [viewer, setViewer] = useState(null);

  useEffect(() => {
    const initializeGridViewer = async () => {
      // Wait for 3Dmol to be loaded
      while (typeof window.$3Dmol === 'undefined') {
        await new Promise(resolve => setTimeout(resolve, 100));
      }

      if (!gridViewerRef.current || viewer) return;

      try {
        console.log(`Creating grid viewer for ${column.viewers.length} molecules`);
        
        // Create a grid viewer for all molecules in this column
        const gridConfig = {
          rows: column.viewers.length,
          cols: 1,
          control_all: false
        };
        
        const gridViewer = window.$3Dmol.createViewerGrid(gridViewerRef.current, gridConfig, {
          backgroundColor: 'transparent',
          antialias: true,
          defaultcolors: window.$3Dmol.rasmolElementColors
        });
        
        // Load each molecule into its grid position
        for (let i = 0; i < column.viewers.length; i++) {
          const molecularData = column.viewers[i];
          if (molecularData.sdfData && molecularData.sdfData.startsWith('file://')) {
            const path = molecularData.sdfData.replace('file://', '');
            try {
              const response = await fetch(`/sdf_files/${path.split('/').pop()}`);
              if (response.ok) {
                const sdfContent = await response.text();
                console.log(`Loading ${molecularData.name} into grid position ${i}`);
                gridViewer.addModel(sdfContent, 'sdf', {}, i);
                gridViewer.setStyle({}, { sphere: { scale: 0.8 } }, i);
                gridViewer.zoomTo(i);
                gridViewer.render(i);
              }
            } catch (err) {
              console.error(`Failed to load ${molecularData.name}:`, err);
            }
          }
        }
        
        setViewer(gridViewer);
      } catch (error) {
        console.error('Failed to create grid viewer:', error);
      }
    };

    initializeGridViewer();
  }, [column.viewers]);

  return (
    <div style={styles.column}>
      <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '20px' }}>
        <h3>{column.query}</h3>
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
      
      {/* Labels for each molecule */}
      {column.viewers.map((mol, idx) => (
        <div key={idx} style={{ 
          background: 'rgba(255, 255, 255, 0.05)', 
          padding: '8px', 
          marginBottom: '5px',
          borderRadius: '4px',
          fontSize: '14px'
        }}>
          {mol.name}
        </div>
      ))}
      
      {/* Grid viewer container */}
      <div 
        ref={gridViewerRef}
        style={{ 
          height: `${column.viewers.length * 200}px`, 
          width: '100%', 
          background: 'transparent',
          border: 'none'
        }}
      />
      
      <div style={{ marginTop: '10px', fontSize: '12px', opacity: 0.7 }}>
        {column.viewers.length} molecule{column.viewers.length !== 1 ? 's' : ''} detected
      </div>
    </div>
  );
};

export default MainLayout;