import React, { useEffect, useRef } from 'react';

const Results = ({ viewers }) => {
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

  return (
    <div className="results-section">
      <div className="snapshots-container">
        {viewers.map((viewer, index) => (
          <MoleculeViewer key={index} viewer={viewer} />
        ))}
      </div>
      <div id="gldiv" ref={glDivRef}></div>
    </div>
  );
};

const MoleculeViewer = ({ viewer }) => {
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
          onClick={() => {
            // TODO: Handle viewer close
          }}
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