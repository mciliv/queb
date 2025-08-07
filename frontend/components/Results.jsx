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
    <div className="molecular-analysis-results-section">
      {/* Molecular analysis progress indicator */}
      {isAnalyzing && (
        <div className="spatial-analysis-header">
          <h3>Analyzing Molecular Content: {getTargetObjectDescription()}</h3>
          <p>Estimating molecules in specified space...</p>
        </div>
      )}
      
      {/* 3dmoljs gridviewer container with main app background */}
      <div id="gldiv" ref={molecularVisualizationContainerRef} className="molecular-gridviewer-container">
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
    <div className="spatial-molecular-column">
      {/* Spatial analysis header */}
      <div className="spatial-target-container">
        <h3 className="target-object-title">{targetObjectName}</h3>
        <div className="detected-molecule-count">
          {molecularViewers.length} molecule{molecularViewers.length !== 1 ? 's' : ''} detected
        </div>
        {molecularViewers.length > 3 && (
          <button 
            className="grid-toggle-btn"
            onClick={() => setUseGridViewer(!useGridViewer)}
            title={useGridViewer ? 'Switch to individual viewers' : 'Switch to grid view'}
          >
            {useGridViewer ? '□' : '⊞'}
          </button>
        )}
      </div>

      {/* 3dmoljs molecular gridviewer with main app background */}
      {useGridViewer ? (
        <div 
          ref={gridViewerRef}
          className="molecule-grid"
          style={{ height: '400px', width: '100%', background: 'transparent' }}
        />
      ) : (
        <div className="threeDmol-molecular-gridviewer">
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
  const [isLoadingMolecularVisualization, setIsLoadingMolecularVisualization] = useState(false);
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

      setIsLoadingMolecularVisualization(true);
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
      } finally {
        setIsLoadingMolecularVisualization(false);
      }
    };

    loadMolecularStructuralData();
  }, [molecularData.sdfData]);

  return (
    <div className="threeDmol-molecular-viewer-container">
      <div className="molecular-viewer-header">
        <span className="detected-molecule-name">{molecularData.name}</span>
        <button 
          className="close-molecular-visualization-btn"
          onClick={onCloseMoleculeVisualization}
          title={`Close ${molecularData.name} visualization (${navigator.platform.toUpperCase().indexOf('MAC') >= 0 ? '⌘' : 'Ctrl'}+W)`}
        >
          ×
        </button>
      </div>
      <div 
        ref={threeDMolViewerRef}
        className="mol-viewer-container"
        style={{ height: '300px', width: '400px', background: 'transparent' }}
      />
      {isLoadingMolecularVisualization && (
        <div className="molecular-visualization-loading">
          <span className="molecular-spinner"></span>
          Loading molecular visualization...
        </div>
      )}
      {molecularVisualizationError && (
        <div className="molecular-visualization-error">
          <div className="molecular-error-message">{molecularVisualizationError}</div>
          {molecularData.smiles && (
            <div className="smiles-structural-fallback">
              SMILES Structure: <code>{molecularData.smiles}</code>
            </div>
          )}
        </div>
      )}
    </div>
  );
};

export default MolecularAnalysisResults;