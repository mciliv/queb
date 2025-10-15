import React, { useState, useEffect, useRef } from 'react';
import logger from '../core/logger.js';

const MoleculeViewer = ({ molecularData }) => {
  const ref = useRef(null);
  const mountRef = useRef(null); // Dedicated imperative mount to avoid React DOM conflicts
  const [status, setStatus] = useState('loading');
  const loggedFailuresRef = useRef(new Set()); // Track logged failures to avoid repetition

  useEffect(() => {
    let cancelled = false;
    const initialize = async () => {
      try {
        // Wait for 3Dmol to load
        while (typeof window.$3Dmol === 'undefined') {
          await new Promise(resolve => setTimeout(resolve, 100));
          if (cancelled) return;
        }
        if (!ref.current) return;
        // Ensure we have a stable child container React does not manage
        if (!mountRef.current) {
          const container = document.createElement('div');
          container.style.width = '100%';
          container.style.height = '100%';
          ref.current.appendChild(container);
          mountRef.current = container;
        }
        const host = mountRef.current;

        // Keep existing viewer visible until new SDF is ready to avoid flash
        const hadViewer = host.childNodes && host.childNodes.length > 0;
        let sdfContent = null;
        if (molecularData.sdfData && molecularData.sdfData.startsWith('file://')) {
          const rawPath = molecularData.sdfData.replace('file://', '');
          // Construct URL properly - rawPath should already be a valid server path like "/sdf_files/filename.sdf"
          const url = rawPath.startsWith('/sdf_files/')
            ? rawPath
            : `/sdf_files/${rawPath}`;
          const response = await fetch(url);
          if (response.ok) {
            sdfContent = await response.text();
            // SDF fetched successfully
          } else {
            const failureKey = `${molecularData.name}-${response.status}`;
            if (!loggedFailuresRef.current.has(failureKey)) {
              logger.info(`SDF fetch failed for ${molecularData.name}: HTTP ${response.status}`);
              logger.debug(`Failed URL: ${url}`);
              loggedFailuresRef.current.add(failureKey);
            }
            // Try fallback approaches for debugging
            if (molecularData.smiles) {
              const sanitizeSmiles = (s) => s.replace(/[^a-zA-Z0-9]/g, ch => ch === '=' ? '__' : '_');
              const fallbackUrl = `/sdf_files/${sanitizeSmiles(molecularData.smiles)}.sdf`;
              const fallbackKey = `${molecularData.name}-fallback`;
              if (!loggedFailuresRef.current.has(fallbackKey)) {
                logger.debug(`Trying fallback URL: ${fallbackUrl}`);
                loggedFailuresRef.current.add(fallbackKey);
              }
              const fallbackResponse = await fetch(fallbackUrl);
              if (fallbackResponse.ok) {
                sdfContent = await fallbackResponse.text();
                logger.info(`SDF fallback successful for ${molecularData.name}`);
              }
            }
          }
        }

        if (sdfContent || molecularData.smiles || molecularData.name) {
          // Now replace the canvas with a fresh viewer
          if (!host) return;
          host.innerHTML = '';
          const viewer = window.$3Dmol.createViewer(host, {
            backgroundColor: '#000000',
            antialias: true,
            defaultcolors: window.$3Dmol.rasmolElementColors
          });
          if (typeof viewer.setBackgroundColor === 'function') {
            viewer.setBackgroundColor(0x000000, 0);
          }
          if (typeof viewer.removeAllModels === 'function') viewer.removeAllModels();
          if (typeof viewer.clear === 'function') viewer.clear();
          if (sdfContent) {
            viewer.addModel(sdfContent, 'sdf');
          } else if (molecularData?.smiles) {
            try { viewer.addModel(molecularData.smiles, 'smiles'); } catch (_) {}
          } else if (molecularData.name) {
            try {
              const encoded = encodeURIComponent(molecularData.name);
              const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encoded}/SDF?record_type=3d`;
              const fetched = await fetch(url);
              if (fetched.ok) {
                const text = await fetched.text();
                viewer.addModel(text, 'sdf');
              }
            } catch (_) {}
          }
          // Use sphere representation with proper element colors
          // For ionic compounds like NaCl, show individual atoms clearly
          viewer.setStyle({}, { 
            sphere: { 
              scale: 0.8
            } 
          });
          // Ensure proper element-based coloring
          viewer.setStyle({element: 'Na'}, { sphere: { color: 'purple', scale: 0.8 } });
          viewer.setStyle({element: 'Cl'}, { sphere: { color: 'green', scale: 0.8 } });
          viewer.setStyle({element: 'H'}, { sphere: { color: 'white', scale: 0.6 } });
          viewer.setStyle({element: 'O'}, { sphere: { color: 'red', scale: 0.8 } });
          viewer.setStyle({element: 'C'}, { sphere: { color: 'gray', scale: 0.8 } });
          viewer.setStyle({element: 'N'}, { sphere: { color: 'blue', scale: 0.8 } });
          if (typeof viewer.resize === 'function') viewer.resize();
          if (typeof viewer.zoomTo === 'function') viewer.zoomTo();
          if (typeof viewer.zoom === 'function') {
            try { viewer.zoom(0.85); } catch (_) {}
          }
          viewer.render();
          setStatus('loaded');
          // Render complete
        } else {
          // If there was already a viewer, keep it and avoid flashing failure UI
          if (!hadViewer) {
            setStatus('failed');
          }
        }
      } catch (error) {
        logger.warn(`Molecule load failed for ${molecularData.name}: ${error.message || 'Unknown error'}`);
        setStatus('failed');
      }
    };

    initialize();
    return () => {
      cancelled = true;
      // Only clear our dedicated mount; avoid touching React-managed nodes
      if (mountRef.current) {
        try { mountRef.current.innerHTML = ''; } catch (_) {}
      }
    };
  }, [molecularData.sdfData, molecularData.name, molecularData.smiles]);

  return (
    <div className="molecule-card">
      <div className="molecule-title">
        {molecularData.name && (
          <a 
            href={`https://en.wikipedia.org/wiki/${encodeURIComponent(molecularData.name)}`}
            target="_blank"
            rel="noopener noreferrer"
          >
            {molecularData.name}
          </a>
        )}
        {status === 'failed' && ' ❌'}
      </div>
      <div 
        ref={ref}
        className="viewer"
      >
        
        {status === 'failed' && '❌'}
      </div>
    </div>
  );
};

export default MoleculeViewer;