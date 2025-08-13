import React, { useState, useEffect, useRef } from 'react';
import { PaymentProvider } from '../components/ui/PaymentContext';
import { useApi } from '../hooks/useApi';
import '../assets/style.css';

// Simple styles object
const styles = {
  app: {
    background: '#000000',
    color: '#ffffff',
    minHeight: '100vh',
    padding: '20px',
    fontFamily: 'system-ui, sans-serif'
  },
  input: {
    width: '100%',
    padding: '12px',
    background: 'rgba(255, 255, 255, 0.1)',
    border: 'none',
    borderRadius: '8px',
    color: '#ffffff',
    fontSize: '14px',
    marginBottom: '10px'
  },
  button: {
    background: 'rgba(255, 255, 255, 0.1)',
    border: 'none',
    borderRadius: '6px',
    color: '#ffffff',
    padding: '8px 12px',
    cursor: 'pointer',
    margin: '5px'
  },
  columns: {
    display: 'flex',
    gap: '20px',
    overflowX: 'auto',
    marginTop: '20px'
  },
  column: {
    minWidth: '400px',
    maxWidth: '400px',
    padding: '20px',
    background: 'rgba(255, 255, 255, 0.02)',
    borderRadius: '8px',
    flexShrink: 0
  },
  header: {
    display: 'flex',
    justifyContent: 'space-between',
    alignItems: 'center',
    marginBottom: '20px',
    padding: '8px 0'
  }
};

// Simple text input component
function TextInput({ value, onChange, onSubmit }) {
  function handleKeyPress(event) {
    if (event.key === 'Enter') {
      onSubmit();
    }
  }

  return (
    <div>
      <input
        style={styles.input}
        type="text"
        placeholder="Enter something to analyze (e.g., water, coffee, salt)"
        value={value}
        onChange={(e) => onChange(e.target.value)}
        onKeyPress={handleKeyPress}
      />
      <button style={styles.button} onClick={onSubmit}>
        Analyze Text
      </button>
    </div>
  );
}

// Simple camera component
function CameraInput({ onImageCapture }) {
  const videoRef = useRef(null);
  const [cameraReady, setCameraReady] = useState(false);

  // Start camera when component loads
  useEffect(() => {
    startCamera();
  }, []);

  async function startCamera() {
    try {
      const stream = await navigator.mediaDevices.getUserMedia({ video: true });
      if (videoRef.current) {
        videoRef.current.srcObject = stream;
        setCameraReady(true);
      }
    } catch (error) {
      console.error('Camera error:', error);
    }
  }

  function captureImage() {
    if (!videoRef.current) return;

    const canvas = document.createElement('canvas');
    const video = videoRef.current;
    canvas.width = video.videoWidth;
    canvas.height = video.videoHeight;
    
    const ctx = canvas.getContext('2d');
    ctx.drawImage(video, 0, 0);
    
    const imageData = canvas.toDataURL('image/jpeg');
    onImageCapture(imageData);
  }

  return (
    <div>
      <video
        ref={videoRef}
        autoPlay
        playsInline
        style={{
          width: '300px',
          height: '200px',
          background: '#333',
          borderRadius: '8px',
          cursor: 'crosshair'
        }}
        onClick={captureImage}
      />
      {cameraReady && (
        <div style={{ marginTop: '10px', fontSize: '12px', opacity: 0.7 }}>
          Click on the camera feed to capture
        </div>
      )}
    </div>
  );
}

// Simple molecule viewer
function MoleculeViewer({ molecule }) {
  const viewerRef = useRef(null);
  const [status, setStatus] = useState('loading');

  useEffect(() => {
    loadMolecule();
  }, [molecule]);

  async function loadMolecule() {
    try {
      // Wait for 3Dmol library to load
      while (!window.$3Dmol) {
        await new Promise(resolve => setTimeout(resolve, 100));
      }

      if (!viewerRef.current) return;

      // Create 3D viewer
      const viewer = window.$3Dmol.createViewer(viewerRef.current, {
        backgroundColor: 'transparent'
      });

      // Load molecule data
      if (molecule.sdfData) {
        const response = await fetch(molecule.sdfData.replace('file://', '/sdf_files/'));
        if (response.ok) {
          const sdfContent = await response.text();
          viewer.addModel(sdfContent, 'sdf');
          viewer.setStyle({}, { sphere: { scale: 0.8 } });
          viewer.zoomTo();
          viewer.render();
          setStatus('loaded');
        }
      }
    } catch (error) {
      console.error('Error loading molecule:', error);
      setStatus('failed');
    }
  }

  return (
    <div style={{ marginBottom: '10px' }}>
      <div style={{ marginBottom: '5px', fontSize: '14px' }}>
        {molecule.name} {status === 'loading' && '⏳'} {status === 'failed' && '❌'}
      </div>
      <div
        ref={viewerRef}
        style={{
          height: '200px',
          background: 'transparent',
          borderRadius: '4px'
        }}
      />
    </div>
  );
}

// Simple column component
function Column({ column, onRemove }) {
  return (
    <div style={styles.column}>
      {/* Column header */}
      <div style={styles.header}>
        <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
          <span>🧪</span>
          <span>{column.query}</span>
          {column.loading && <span>⏳</span>}
        </div>
        <button style={styles.button} onClick={onRemove}>
          ×
        </button>
      </div>

      {/* Show error if analysis failed */}
      {column.failed && (
        <div style={{
          background: 'rgba(255, 193, 7, 0.1)',
          border: '1px solid rgba(255, 193, 7, 0.3)',
          padding: '10px',
          borderRadius: '4px',
          marginBottom: '10px',
          fontSize: '13px'
        }}>
          💡 AI failed to analyze - please report this
        </div>
      )}

      {/* Show loading message */}
      {column.loading && (!column.viewers || column.viewers.length === 0) && (
        <div style={{ textAlign: 'center', padding: '20px', opacity: 0.5 }}>
          Analyzing molecules...
        </div>
      )}

      {/* Show molecules */}
      {column.viewers && column.viewers.map((molecule, index) => (
        <MoleculeViewer key={index} molecule={molecule} />
      ))}
    </div>
  );
}

// Main app component
function App() {
  // State variables (data that can change)
  const [textInput, setTextInput] = useState('');
  const [showCamera, setShowCamera] = useState(false);
  const [columns, setColumns] = useState([]);
  const [isAnalyzing, setIsAnalyzing] = useState(false);

  // API functions
  const { analyzeText, analyzeImage, generateSDFs } = useApi();

  // Load 3D molecule library
  useEffect(() => {
    if (!document.getElementById('3dmol-script')) {
      const script = document.createElement('script');
      script.id = '3dmol-script';
      script.src = 'https://3Dmol.org/build/3Dmol-min.js';
      document.head.appendChild(script);
    }
  }, []);

  // Analyze text input
  async function handleTextAnalysis() {
    if (!textInput.trim() || isAnalyzing) return;

    setIsAnalyzing(true);

    // Create new column
    const newColumn = {
      id: Date.now(),
      query: textInput,
      viewers: [],
      loading: true,
      failed: false
    };
    setColumns(prev => [...prev, newColumn]);

    try {
      // Get analysis from AI
      const result = await analyzeText(textInput);
      const molecules = result.molecules || [];

      if (molecules.length > 0) {
        // Generate 3D structure files
        const smiles = molecules.map(mol => mol.smiles).filter(Boolean);
        const sdfResult = await generateSDFs(smiles, false);

        // Create molecule objects with 3D data
        const viewerData = molecules.map((mol, index) => ({
          name: mol.name || textInput,
          smiles: mol.smiles,
          sdfData: sdfResult.sdfPaths?.[index] ? `file://${sdfResult.sdfPaths[index]}` : null
        }));

        // Update column with results
        setColumns(prev => prev.map(col =>
          col.id === newColumn.id
            ? { ...col, viewers: viewerData, loading: false, failed: false }
            : col
        ));
      } else {
        // No molecules found
        setColumns(prev => prev.map(col =>
          col.id === newColumn.id
            ? { ...col, loading: false, failed: true }
            : col
        ));
      }
    } catch (error) {
      console.error('Analysis failed:', error);
      // Mark column as failed
      setColumns(prev => prev.map(col =>
        col.id === newColumn.id
          ? { ...col, loading: false, failed: true }
          : col
      ));
    }

    setIsAnalyzing(false);
    setTextInput(''); // Clear input after analysis
  }

  // Analyze image from camera
  async function handleImageAnalysis(imageData) {
    setIsAnalyzing(true);

    // Create new column
    const newColumn = {
      id: Date.now(),
      query: 'Camera capture',
      viewers: [],
      loading: true,
      failed: false
    };
    setColumns(prev => [...prev, newColumn]);

    try {
      // Get analysis from AI
      const result = await analyzeImage(imageData, 'Camera capture');
      const molecules = result.molecules || [];

      if (molecules.length > 0) {
        // Generate 3D structure files
        const smiles = molecules.map(mol => mol.smiles).filter(Boolean);
        const sdfResult = await generateSDFs(smiles, false);

        // Create molecule objects with 3D data
        const viewerData = molecules.map((mol, index) => ({
          name: mol.name || 'Camera capture',
          smiles: mol.smiles,
          sdfData: sdfResult.sdfPaths?.[index] ? `file://${sdfResult.sdfPaths[index]}` : null
        }));

        // Update column with results
        setColumns(prev => prev.map(col =>
          col.id === newColumn.id
            ? { ...col, viewers: viewerData, loading: false, failed: false }
            : col
        ));
      } else {
        // No molecules found
        setColumns(prev => prev.map(col =>
          col.id === newColumn.id
            ? { ...col, loading: false, failed: true }
            : col
        ));
      }
    } catch (error) {
      console.error('Image analysis failed:', error);
      // Mark column as failed
      setColumns(prev => prev.map(col =>
        col.id === newColumn.id
          ? { ...col, loading: false, failed: true }
          : col
      ));
    }

    setIsAnalyzing(false);
  }

  // Remove a column
  function removeColumn(columnId) {
    setColumns(prev => prev.filter(col => col.id !== columnId));
  }

  return (
    <PaymentProvider config={{ enabled: false, devMode: true }}>
      <div style={styles.app}>
        <h1 style={{ marginBottom: '30px' }}>🧪 Structuralizer</h1>

        {/* Text input section */}
        <div style={{ marginBottom: '20px', maxWidth: '600px' }}>
          <TextInput
            value={textInput}
            onChange={setTextInput}
            onSubmit={handleTextAnalysis}
          />
        </div>

        {/* Camera toggle */}
        <div style={{ marginBottom: '20px' }}>
          <button
            style={styles.button}
            onClick={() => setShowCamera(!showCamera)}
          >
            {showCamera ? 'Hide Camera' : 'Show Camera'}
          </button>
        </div>

        {/* Camera section */}
        {showCamera && (
          <div style={{ marginBottom: '20px' }}>
            <CameraInput onImageCapture={handleImageAnalysis} />
          </div>
        )}

        {/* Results columns */}
        <div style={styles.columns}>
          {columns.map(column => (
            <Column
              key={column.id}
              column={column}
              onRemove={() => removeColumn(column.id)}
            />
          ))}
        </div>
      </div>
    </PaymentProvider>
  );
}

export default App;