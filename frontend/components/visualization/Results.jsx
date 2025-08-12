import React, { useEffect, useRef, useState } from 'react';

const MolecularAnalysisResults = ({ moleculeViewers, setMoleculeViewers, lastAnalysisResult, isAnalyzing, currentAnalysisType, targetObjectInput }) => {
  const molecularVisualizationContainerRef = useRef(null);

  useEffect(() => {
    // Load 3Dmol.js library for molecular visualization
    if (typeof window.$3Dmol === 'undefined') {
      const threeDMolScript = document.createElement('script');
      threeDMolScript.src = 'https://3Dmol.org/build/3Dmol-min.js';
      threeDMolScript.async = true;
      threeDMolScript.onload = () => {
        console.log('3Dmol.js library loaded successfully');
      };
      threeDMolScript.onerror = () => {
        console.error('Failed to load 3Dmol.js library');
      };
      document.head.appendChild(threeDMolScript);
    } else {
      console.log('3Dmol.js already loaded');
    }
  }, []);

  // Handle keyboard shortcuts for closing molecular viewers
  useEffect(() => {
    const handleMolecularViewerKeyboard = (event) => {
      // Ignore shortcuts when user is typing in input fields
      if (event.target.tagName === 'INPUT' || event.target.tagName === 'TEXTAREA') {
        return;
      }

      const isMacOS = navigator.platform.toUpperCase().indexOf('MAC') >= 0;
      const modifierKey = isMacOS ? event.metaKey : event.ctrlKey;

      if (modifierKey && event.key === 'w' && moleculeViewers.length > 0) {
        event.preventDefault();
        // Close the most recent molecular viewer
        setMoleculeViewers(previousViewers => previousViewers.slice(0, -1));
      }
    };

    document.addEventListener('keydown', handleMolecularViewerKeyboard);
    return () => document.removeEventListener('keydown', handleMolecularViewerKeyboard);
  }, [moleculeViewers.length, setMoleculeViewers]);

  // Extract target object description from analysis results
  const getTargetObjectDescription = () => {
    if (isAnalyzing) {
      // During molecular analysis, show the input target object if available
      return targetObjectInput || 'Unknown Target';
    }
    
    if (lastAnalysisResult) {
      // For completed analysis, get the target object description from the result
      // Check multiple possible locations for the object description
      if (lastAnalysisResult.object) {
        return lastAnalysisResult.object;
      }
      if (lastAnalysisResult.result && lastAnalysisResult.result.object) {
        return lastAnalysisResult.result.object;
      }
      if (lastAnalysisResult.input) {
        return lastAnalysisResult.input;
      }
      // For image/camera analysis, the result might be directly the analysis result
      if (lastAnalysisResult.molecules && lastAnalysisResult.molecules.length > 0) {
        // This is likely an image analysis result
        return lastAnalysisResult.object || 'Captured Object';
      }
    }
    
    return 'Molecular Analysis';
  };

  // Group molecular viewers by analysis session (assuming they come from the same analysis if added together)
  const molecularAnalysisGroups = [];
  if (lastAnalysisResult && moleculeViewers.length > 0) {
    // For now, treat all current viewers as one analysis group
    // In the future, could enhance to track multiple analysis sessions
    molecularAnalysisGroups.push({
      targetObjectName: getTargetObjectDescription(),
      molecularViewers: moleculeViewers,
      analysisResult: lastAnalysisResult
    });
  }

  return (
    <div style={{ width: '100%', background: '#000' }}>
      <div
        id="gldiv"
        ref={molecularVisualizationContainerRef}
        style={{
          display: 'flex',
          flexDirection: 'row',
          gap: '20px',
          padding: '20px',
          background: '#000',
          minHeight: '100px',
          width: '100%',
          overflowX: 'auto',
          overflowY: 'hidden'
        }}
      >
        {molecularAnalysisGroups.map((spatialAnalysisGroup, groupIndex) => (
          <SpatialMolecularColumn
            key={groupIndex}
            targetObjectName={spatialAnalysisGroup.targetObjectName}
            molecularViewers={spatialAnalysisGroup.molecularViewers}
            spatialAnalysisResult={spatialAnalysisGroup.analysisResult}
            onCloseMolecularViewer={(viewerIndex) => {
              setMoleculeViewers(previousViewers => previousViewers.filter((_, i) => i !== viewerIndex));
            }}
          />
        ))}
      </div>
    </div>
  );
};

const SpatialMolecularColumn = ({ targetObjectName, molecularViewers, spatialAnalysisResult, onCloseMolecularViewer }) => {
  const gridViewerRef = useRef(null);
  const [useGridViewer, setUseGridViewer] = useState(false);

  useEffect(() => {
    // Use 3Dmol.js createViewerGrid for multiple molecules as per documentation
    if (useGridViewer && molecularViewers.length > 1 && window.$3Dmol && gridViewerRef.current) {
      const gridConfig = {
        rows: Math.ceil(Math.sqrt(molecularViewers.length)),
        cols: Math.ceil(Math.sqrt(molecularViewers.length)),
        control_all: true  // Synchronize all viewers in grid
      };
      
      try {
        const gridViewer = window.$3Dmol.createViewerGrid(gridViewerRef.current, gridConfig, {
          backgroundColor: 'transparent',  // Transparent background as in original
          antialias: true,
          defaultcolors: window.$3Dmol.rasmolElementColors,
          ui: false,
          showControls: false,
          showInfo: false
        });
        
        // Add each molecule to the grid
        molecularViewers.forEach((molecularData, index) => {
          if (molecularData.sdfData) {
            gridViewer.addModel(molecularData.sdfData, 'sdf', {}, index);
            // Use original sphere-only styling with proper van der Waals radii
            gridViewer.setStyle({}, {
              sphere: { scale: 0.8 }  // Van der Waals radii as in original
            }, index);
          }
        });
        
        gridViewer.zoomTo();
        gridViewer.render();
      } catch (error) {
        console.error('Failed to create 3Dmol.js grid viewer:', error);
        // Don't set state in useEffect to avoid loops - let user manually toggle
      }
    }
  }, [useGridViewer, molecularViewers]);

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: '20px', minWidth: '400px', flexShrink: 0 }}>
      <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
        {molecularViewers.length > 3 && (
          <button 
            onClick={() => setUseGridViewer(!useGridViewer)}
            style={{
              background: 'rgba(255,255,255,0.1)',
              border: 'none',
              color: 'rgba(255,255,255,0.8)',
              fontSize: '16px',
              padding: '8px',
              borderRadius: '4px',
              cursor: 'pointer'
            }}
            aria-label="Toggle grid view"
          >
            {useGridViewer ? '□' : '⊞'}
          </button>
        )}
      </div>

      {useGridViewer ? (
        <div 
          ref={gridViewerRef}
          style={{ height: '400px', width: '100%', background: 'transparent' }}
        />
      ) : (
        <div style={{ display: 'flex', flexDirection: 'column', gap: '20px', width: '100%' }}>
          {molecularViewers.map((molecularViewer, viewerIndex) => (
            <ThreeDimensionalMolecularViewer 
              key={viewerIndex} 
              molecularData={molecularViewer} 
              spatialIndex={viewerIndex}
              onCloseMoleculeVisualization={() => onCloseMolecularViewer(viewerIndex)}
              isWithinGridLayout={true}
            />
          ))}
        </div>
      )}
    </div>
  );
};

const ThreeDimensionalMolecularViewer = ({ molecularData, spatialIndex, onCloseMoleculeVisualization, isWithinGridLayout = false }) => {
  const threeDMolViewerRef = useRef(null);
  const [structuralDataFile, setStructuralDataFile] = useState(null);
  const [molecularVisualizationError, setMolecularVisualizationError] = useState(null);

  useEffect(() => {
    const loadMolecularStructuralData = async () => {
      if (!molecularData.sdfData || !threeDMolViewerRef.current) {
        return;
      }

      // Wait for 3Dmol.js to load if not available yet
      if (!window.$3Dmol) {
        console.log('Waiting for 3Dmol.js to load...');
        const waitFor3Dmol = () => {
          return new Promise((resolve) => {
            const checkInterval = setInterval(() => {
              if (window.$3Dmol) {
                clearInterval(checkInterval);
                console.log('3Dmol.js is now available');
                resolve();
              }
            }, 100);
            // Timeout after 10 seconds
            setTimeout(() => {
              clearInterval(checkInterval);
              console.error('Timeout waiting for 3Dmol.js to load');
              resolve();
            }, 10000);
          });
        };
        await waitFor3Dmol();
        
        if (!window.$3Dmol) {
          setMolecularVisualizationError('3Dmol.js library failed to load');
          return;
        }
      }

      setMolecularVisualizationError(null);

      try {
        let sdfStructuralContent = molecularData.sdfData;

        // If sdfData is a file path, fetch the structural content
        if (typeof sdfStructuralContent === 'string' && sdfStructuralContent.startsWith('file://')) {
          const structuralFilePath = sdfStructuralContent.replace('file://', '');
          const structuralDataResponse = await fetch(`/sdf_files/${structuralFilePath.split('/').pop()}`);
          if (!structuralDataResponse.ok) {
            throw new Error('Failed to load molecular structural data file');
          }
          sdfStructuralContent = await structuralDataResponse.text();
        }

        // If we have SDF structural content, render 3D molecular visualization using 3Dmol.js
        if (sdfStructuralContent && sdfStructuralContent.trim()) {
          console.log('Creating 3Dmol viewer with SDF data length:', sdfStructuralContent.length);
          console.log('Container element:', threeDMolViewerRef.current);
          console.log('3Dmol available:', !!window.$3Dmol);
          
          const threeDMolViewer = window.$3Dmol.createViewer(threeDMolViewerRef.current, {
            backgroundColor: 'transparent',  // Transparent background as in original
            antialias: true,
            defaultcolors: window.$3Dmol.rasmolElementColors,
            // Disable 3Dmol.js UI to match app's "stupid simple" design
            ui: false,
            showControls: false,
            showInfo: false
          });
          
          console.log('3Dmol viewer created:', threeDMolViewer);
          
          // Add molecular model with SDF format support
          const model = threeDMolViewer.addModel(sdfStructuralContent, 'sdf');
          console.log('Model added:', model);
          
          // CRITICAL: Use ONLY sphere representation with van der Waals radii at 0.8 scale (as in original)
          threeDMolViewer.setStyle({}, { 
            sphere: { 
              scale: 0.8  // Van der Waals radii scaling as in original implementation
            } 
          });
          
          // Optimize view and render (as in original)
          threeDMolViewer.zoomTo();
          threeDMolViewer.render();
          
          console.log('3Dmol viewer rendered successfully');
          setStructuralDataFile(sdfStructuralContent);
        } else {
          console.error('No SDF structural data available for visualization');
          setMolecularVisualizationError('No molecular structural data available');
        }
      } catch (molecularRenderingError) {
        console.error('Failed to load molecular structural data:', molecularRenderingError);
        setMolecularVisualizationError('Failed to load molecule visualization');
      }
    };

    loadMolecularStructuralData();
  }, [molecularData.sdfData]);

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: '8px', background: 'rgba(255,255,255,0.02)', borderRadius: '4px', padding: '12px' }}>
      <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', gap: '12px', marginBottom: '8px' }}>
        <span style={{ color: 'rgba(255,255,255,0.9)', fontSize: '14px', fontWeight: 500, flex: 1, textAlign: 'left' }}>
          {molecularData.name}
        </span>
        <button 
          onClick={onCloseMoleculeVisualization}
          style={{
            background: 'none',
            border: 'none',
            color: 'rgba(255,255,255,0.6)',
            fontSize: '18px',
            cursor: 'pointer',
            padding: '4px 8px',
            borderRadius: '4px',
            flexShrink: 0
          }}
          aria-label={`Close ${molecularData.name}`}
        >
          ×
        </button>
      </div>
      <div 
        ref={threeDMolViewerRef}
        style={{
          position: 'relative',
          width: '400px',
          height: '300px',
          background: 'transparent',
          borderRadius: '4px',
          overflow: 'hidden',
          margin: 0,
          border: 'none',
          flexShrink: 0,
          display: 'block'
        }}
      />

      {molecularVisualizationError && (
        <div style={{ padding: '16px', background: 'rgba(255, 50, 50, 0.1)', borderRadius: '4px' }}>
          <div style={{ color: 'rgba(255, 255, 255, 0.8)', fontSize: '13px', marginBottom: '8px' }}>{molecularVisualizationError}</div>
          {molecularData.smiles && (
            <div style={{ color: 'rgba(255, 255, 255, 0.5)', fontSize: '11px', fontFamily: 'monospace', background: 'rgba(255,255,255,0.05)', padding: '8px', borderRadius: '3px', wordBreak: 'break-all' }}>
              SMILES Structure: <code style={{ color: 'rgba(255,255,255,0.7)', background: 'none' }}>{molecularData.smiles}</code>
            </div>
          )}
        </div>
      )}
    </div>
  );
};

export default MolecularAnalysisResults;