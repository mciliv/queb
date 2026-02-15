/**
 * Hook Migration Example - How to refactor from useApi to focused hooks
 * 
 * This example shows how to migrate from the generic useApi hook
 * to the new user-focused compositional hooks.
 */

// ============================================
// BEFORE: Using generic useApi
// ============================================

import { useApi } from '../hooks/useApi';

function OldTextAnalysisComponent() {
  const [text, setText] = useState('');
  const [results, setResults] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  
  const { structuresFromText, generateSDFs } = useApi();
  
  const handleAnalyze = async () => {
    setLoading(true);
    setError(null);
    
    try {
      // Generic API call with unclear intent
      const result = await structuresFromText(text, 'GPT-5');
      const molecules = result.molecules || result.chemicals || [];
      
      // Need to manually generate structures
      if (molecules.length > 0) {
        const smilesArray = molecules.map(m => m.smiles).filter(Boolean);
        const sdfResult = await generateSDFs(smilesArray);
        // Complex mapping logic...
      }
      
      setResults(molecules);
    } catch (err) {
      setError(err.message);
    } finally {
      setLoading(false);
    }
  };
  
  return (
    <div>
      <input value={text} onChange={e => setText(e.target.value)} />
      <button onClick={handleAnalyze} disabled={loading}>
        Analyze
      </button>
      {/* Display results */}
    </div>
  );
}

// ============================================
// AFTER: Using focused hooks
// ============================================

import { useAnalyzeText, useMolecule3D } from '../hooks';

function NewTextAnalysisComponent() {
  const [text, setText] = useState('');
  
  // Clear intent: "I want to analyze text"
  const { 
    analyze, 
    molecules, 
    isAnalyzing, 
    error 
  } = useAnalyzeText();
  
  // Clear intent: "I want 3D structures"
  const { generateMultiple3D } = useMolecule3D();
  
  const handleAnalyze = async () => {
    // Much simpler - the hook handles all complexity
    const result = await analyze(text, { mode: 'database' });
    
    // Generate 3D structures if needed
    if (result.molecules.length > 0) {
      await generateMultiple3D(result.molecules);
    }
  };
  
  return (
    <div>
      <input value={text} onChange={e => setText(e.target.value)} />
      <button onClick={handleAnalyze} disabled={isAnalyzing}>
        Analyze
      </button>
      
      {/* State is managed by the hook */}
      {error && <div className="error">{error}</div>}
      
      {molecules.map(mol => (
        <div key={mol.name}>
          {mol.name} - {mol.formula}
        </div>
      ))}
    </div>
  );
}

// ============================================
// CAMERA EXAMPLE
// ============================================

// BEFORE: Complex camera handling
function OldCameraComponent() {
  const videoRef = useRef();
  const [stream, setStream] = useState(null);
  const [processing, setProcessing] = useState(false);
  const { structuralizeImage } = useApi();
  
  const startCamera = async () => {
    try {
      const stream = await navigator.mediaDevices.getUserMedia({/*...*/});
      videoRef.current.srcObject = stream;
      setStream(stream);
    } catch (err) {
      // Handle error
    }
  };
  
  const handleClick = async (e) => {
    // Calculate coordinates
    // Capture frame
    // Convert to base64
    // Send to API
    const result = await structuralizeImage(base64, x, y, /*...*/);
    // Process result...
  };
  
  // Cleanup code...
}

// AFTER: Simple camera hook
function NewCameraComponent() {
  const videoRef = useRef();
  const { 
    startCamera, 
    stopCamera, 
    analyzeClick,
    molecules,
    isAnalyzing,
    identifiedObject
  } = useAnalyzeCamera();
  
  useEffect(() => {
    // Hook handles all camera setup
    startCamera(videoRef.current);
    
    return () => stopCamera();
  }, []);
  
  const handleClick = async (e) => {
    const rect = videoRef.current.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;
    
    // One simple call - hook handles everything
    await analyzeClick(x, y);
  };
  
  return (
    <div>
      <video ref={videoRef} onClick={handleClick} />
      {isAnalyzing && <div>Analyzing...</div>}
      {identifiedObject && <h3>Found: {identifiedObject}</h3>}
      {molecules.map(mol => <MoleculeViewer key={mol.name} molecule={mol} />)}
    </div>
  );
}

// ============================================
// KEY BENEFITS
// ============================================

/**
 * 1. CLEARER INTENT
 *    - useAnalyzeText → "I want to analyze text"
 *    - useAnalyzeCamera → "I want to analyze from camera"
 *    - No generic "api.doSomething" calls
 * 
 * 2. BUILT-IN STATE MANAGEMENT
 *    - Hooks manage loading, error, and result states
 *    - No manual state tracking needed
 * 
 * 3. SIMPLIFIED ERROR HANDLING
 *    - Errors are caught and exposed by hooks
 *    - Consistent error interface across all hooks
 * 
 * 4. FOCUSED FUNCTIONALITY
 *    - Each hook does one thing well
 *    - No mixing of concerns (analysis vs 3D generation)
 * 
 * 5. BETTER TYPESCRIPT SUPPORT
 *    - Each hook can have specific types
 *    - No generic "any" API responses
 */


