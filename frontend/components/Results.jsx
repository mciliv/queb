import React, { useEffect, useRef } from 'react';

const Results = ({ viewers, setViewers }) => {
  const glDivRef = useRef(null);

  useEffect(() => {
    // Load 3Dmol.js dynamically if not already loaded
    if (typeof window.$3Dmol === 'undefined') {
      const script = document.createElement('script');
      script.src = 'https://3Dmol.org/build/3Dmol-min.js';
      script.async = true;
      document.head.appendChild(script);
    }
  }, []);

  // Handle keyboard shortcuts for closing viewers
  useEffect(() => {
    const handleKeyDown = (event) => {
      // Ignore shortcuts when typing in input fields
      if (event.target.tagName === 'INPUT' || event.target.tagName === 'TEXTAREA') {
        return;
      }

      const isMac = navigator.platform.toUpperCase().indexOf('MAC') >= 0;
      const modifier = isMac ? event.metaKey : event.ctrlKey;

      if (modifier && event.key === 'w' && viewers.length > 0) {
        event.preventDefault();
        // Close the last viewer
        setViewers(prev => prev.slice(0, -1));
      }
    };

    document.addEventListener('keydown', handleKeyDown);
    return () => document.removeEventListener('keydown', handleKeyDown);
  }, [viewers.length, setViewers]);

  return (
    <div className="results-section">
      <div className="snapshots-container">
        {viewers.map((viewer, index) => (
          <MoleculeViewer 
            key={index} 
            viewer={viewer} 
            index={index}
            onClose={() => setViewers(prev => prev.filter((_, i) => i !== index))}
          />
        ))}
      </div>
      <div id="gldiv" ref={glDivRef}></div>
    </div>
  );
};

const MoleculeViewer = ({ viewer, index, onClose }) => {
  const viewerRef = useRef(null);

  useEffect(() => {
    if (viewer.sdfData && viewerRef.current && window.$3Dmol) {
      const viewer3d = window.$3Dmol.createViewer(viewerRef.current, {
        backgroundColor: '#1a1a1a'
      });
      
      viewer3d.addModel(viewer.sdfData, 'sdf');
      viewer3d.setStyle({}, {
        sphere: { 
          scale: 0.8,
          colorscheme: 'Jmol'
        }
      });
      viewer3d.zoomTo();
      viewer3d.render();
    }
  }, [viewer.sdfData]);

  return (
    <div className="molecule-viewer">
      <div className="viewer-header">
        <span className="molecule-name">{viewer.name}</span>
        <button 
          className="close-btn"
          onClick={onClose}
          title={`Close ${viewer.name} (${navigator.platform.toUpperCase().indexOf('MAC') >= 0 ? '⌘' : 'Ctrl'}+W)`}
        >
          ×
        </button>
      </div>
      <div 
        ref={viewerRef}
        className="viewer-container"
        style={{ height: '300px', width: '100%' }}
      />
    </div>
  );
};

export default Results;