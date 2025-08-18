import React, { useState, useEffect, useCallback, useRef } from 'react';
import { PaymentProvider, usePayment } from '../components/ui/PaymentContext';
import { useApi } from '../hooks/useApi';
import { PRESET_VISUAL_TESTS, TEST_MOLECULES, SMILES_NAME_MAP } from './constants.js';
import '../assets/style.css';

const isMobileDevice = () => {
  return /Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini/i.test(navigator.userAgent) ||
         (window.innerWidth <= 768 && 'ontouchstart' in window);
};

// Configuration
const PAYMENT_CONFIG = {
  enabled: false, // Set to true to enable payment functionality
  devMode: window.location.hostname === 'localhost' || window.location.hostname === '127.0.0.1'
};

// Styles moved to CSS in assets/style.css
const styles = {
  column: {
    minWidth: '400px',
    background: 'transparent',
    padding: '20px',
    position: 'relative',
    zIndex: 1,
    userSelect: 'text'
  }
};

// Input Components (Top Row)
const TextInput = ({ value, onChange, onSubmit, isProcessing, error }) => {
  const [localError, setLocalError] = useState('');
  const [isValidating, setIsValidating] = useState(false);
  const inputRef = useRef(null);
  const lastTriggerTimeRef = useRef(0);
  const isMac = navigator.userAgent.toLowerCase().includes('mac');
  const keyboardHint = isMac ? '⌘K' : 'Ctrl+K';

  useEffect(() => {
    if (localError && value) {
      setLocalError('');
    }
  }, [value, localError]);

  useEffect(() => {
    if (inputRef.current) {
      inputRef.current.focus();
    }
  }, []);

  const validateInput = (text) => {
    if (!text || !text.trim()) {
      return 'Enter a thing to structuralize';
    }
    
    const trimmed = text.trim().toLowerCase();
    
    if (trimmed.length < 2) {
      return 'Input must be at least 2 characters long';
    }
    if (trimmed.length > 500) {
      return 'Input must be less than 500 characters';
    }
    
    const nonPhysicalPatterns = [
      /^(love|hate|happy|sad|angry|joy|fear|hope|dream|idea|thought|feeling|emotion)/,
      /^(running|walking|talking|thinking|sleeping|eating|drinking)$/,
      /^[a-z]{1,2}$/,
      /^[^a-z]*$/,
      /^(asdf|qwerty|test|random|nothing|something|anything|everything)$/i,
      /^(the|a|an|and|or|but|if|then|when|where|why|how|what|who)$/
    ];
    
    if (nonPhysicalPatterns.some(pattern => pattern.test(trimmed))) {
      return 'Please describe a real, physical object (food, materials, plants, etc.)';
    }
    
    return null;
  };

  const handleSubmit = async () => {
    const validationError = validateInput(value);
    if (validationError) {
      setLocalError(validationError);
      return;
    }

    setIsValidating(true);
    try {
      await onSubmit(value.trim());
      setLocalError('');
    } catch (err) {
      setLocalError(err.message || 'Structuralization failed. Please try again.');
    } finally {
      setIsValidating(false);
    }
  };

  const handleKeyDown = async (e) => {
    if (e.key === 'Enter') {
      e.preventDefault();
      e.stopPropagation();
      await handleSubmit();
    }
  };

  const displayError = localError || error;

  return (
    <div className="input-wrapper">
      <div className="input-row">
        <input
          ref={inputRef}
          id="object-input"
          type="text"
          placeholder="Type to analyze molecules..."
          className={`input-base${displayError ? ' input-error' : ''}`}
          value={value}
          onChange={(e) => onChange(e.target.value)}
          onKeyDown={handleKeyDown}
          aria-describedby={displayError ? 'input-error' : undefined}
        />
        
        {value.trim() && (
          <button 
            className="btn-icon"
            onClick={handleSubmit}
            aria-label="Structuralize"
          >
            →
         </button>
        )}
        
        {/* Modern keyboard hint badge */}
        {!value.trim() && (
          <div className="kbd-hint">{keyboardHint}</div>
        )}
      </div>
      {displayError && (
        <div id="input-error" className="error-text" role="alert">
          {displayError}
        </div>
      )}
    </div>
  );
};

const ModeSelector = ({ cameraMode, setCameraMode, photoMode, setPhotoMode, linkMode, setLinkMode }) => {
  const handleModeSelect = (mode) => {
    setCameraMode(false);
    setPhotoMode(false);
    if (setLinkMode) setLinkMode(false);
    
    switch(mode) {
      case 'camera':
        setCameraMode(true);
        break;
      case 'photo':
        setPhotoMode(true);
        break;
      case 'link':
        if (setLinkMode) setLinkMode(true);
        break;
    }
  };

  const buttonStyle = {
    background: 'transparent',
    border: 'none',
    color: '#ffffff',
    padding: '8px 12px',
    borderRadius: '4px',
    cursor: 'pointer',
    display: 'flex',
    alignItems: 'center',
    fontSize: '13px',
    transition: 'all 0.2s'
  };

  const activeStyle = {
    background: 'rgba(255, 255, 255, 0.1)'
  };

  return (
    <div className="mode-row">
      <button 
        style={{...buttonStyle, ...(cameraMode ? activeStyle : {})}}
        onClick={() => handleModeSelect('camera')}
        title="Capture from camera"
      >
        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
          <circle cx="12" cy="12" r="5" fill="currentColor" opacity="0.8"/>
          <circle cx="12" cy="12" r="9" stroke="currentColor" fill="none"/>
        </svg>
      </button>

      <button 
        style={{...buttonStyle, ...(photoMode ? activeStyle : {})}}
        onClick={() => handleModeSelect('photo')}
        title="Upload image"
      >
        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
          <rect x="3" y="3" width="18" height="18" rx="2" ry="2"/>
          <circle cx="9" cy="9" r="2"/>
          <path d="M21 15l-3.086-3.086a2 2 0 0 0-2.828 0L6 21"/>
        </svg>
      </button>

      <button 
        style={{...buttonStyle, ...(linkMode ? activeStyle : {})}}
        onClick={() => handleModeSelect('link')}
        title="Enter image link"
      >
        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
          <path d="M10 13a5 5 0 0 0 7.54.54l3-3a5 5 0 0 0-7.07-7.07l-1.72 1.71"/>
          <path d="M14 11a5 5 0 0 0-7.54-.54l-3 3a5 5 0 0 0 7.07 7.07l1.71-1.71"/>
        </svg>
      </button>
    </div>
  );
};

const CameraSection = ({ isProcessing, setIsProcessing, setCurrentAnalysisType, onAnalysisComplete }) => {
  const videoRef = useRef(null);
  const [hasPermission, setHasPermission] = useState(false);
  const [permissionMessage, setPermissionMessage] = useState('');
  const [showSwitchCamera, setShowSwitchCamera] = useState(false);
  const [stream, setStream] = useState(null);
  const [clickPosition, setClickPosition] = useState(null);
  const [showOutline, setShowOutline] = useState(false);
  const { checkPaymentRequired } = usePayment();
  const { analyzeImage } = useApi();

  useEffect(() => {
    requestCameraAccess();
    return () => {
      if (stream) {
        stream.getTracks().forEach(track => track.stop());
      }
    };
  }, []);

  const requestCameraAccess = async () => {
    try {
      const mediaStream = await navigator.mediaDevices.getUserMedia({ 
        video: { facingMode: 'environment' },
        audio: false 
      });
      
      if (videoRef.current) {
        videoRef.current.srcObject = mediaStream;
      }
      
      setStream(mediaStream);
      setHasPermission(true);
      
      const devices = await navigator.mediaDevices.enumerateDevices();
      const videoDevices = devices.filter(device => device.kind === 'videoinput');
      setShowSwitchCamera(videoDevices.length > 1);
    } catch (error) {
      console.error('Camera access error:', error);
      setPermissionMessage('Camera access required to structuralize from camera');
      setHasPermission(false);
    }
  };

  const handleCameraClick = async (event) => {
    if (!hasPermission) {
      await requestCameraAccess();
      return;
    }
    if (!videoRef.current || isProcessing) return;

    if (checkPaymentRequired()) {
      return;
    }

    // Get click position relative to video element
    const rect = event.currentTarget.getBoundingClientRect();
    const x = event.clientX - rect.left;
    const y = event.clientY - rect.top;
    
    // Store click position for outline box
    setClickPosition({ x, y });
    setShowOutline(true);
    
    // Hide outline after 2 seconds
    setTimeout(() => {
      setShowOutline(false);
    }, 2000);

    setIsProcessing(true);
    setCurrentAnalysisType('camera');
    
    try {
      const canvas = document.createElement('canvas');
      const video = videoRef.current;
      const width = video.videoWidth;
      const height = video.videoHeight;
      canvas.width = width;
      canvas.height = height;
      const ctx = canvas.getContext('2d');

      // If user clicked, draw a centered crop square around click (25% of min dimension)
      const cropSide = Math.floor(Math.min(width, height) * 0.25);
      const cropX = Math.max(0, Math.min(width - cropSide, Math.round((x * (width / rect.width)) - cropSide / 2)));
      const cropY = Math.max(0, Math.min(height - cropSide, Math.round((y * (height / rect.height)) - cropSide / 2)));

      ctx.drawImage(video, 0, 0);
      const fullImageDataUrl = canvas.toDataURL('image/jpeg', 0.8);

      // Extract cropped region as separate base64 to boost accuracy
      const cropCanvas = document.createElement('canvas');
      cropCanvas.width = cropSide;
      cropCanvas.height = cropSide;
      const cropCtx = cropCanvas.getContext('2d');
      cropCtx.drawImage(canvas, cropX, cropY, cropSide, cropSide, 0, 0, cropSide, cropSide);
      const croppedImageDataUrl = cropCanvas.toDataURL('image/jpeg', 0.8);
      
      // Pass click coordinates to the API (scaled to video dimensions)
      const scaleX = video.videoWidth / rect.width;
      const scaleY = video.videoHeight / rect.height;
      const centerX = Math.round(x * scaleX);
      const centerY = Math.round(y * scaleY);
      
      const result = await analyzeImage(fullImageDataUrl, 'Camera capture', centerX, centerY, cropX + Math.floor(cropSide/2), cropY + Math.floor(cropSide/2), cropSide);
      
      if (onAnalysisComplete) {
        onAnalysisComplete(result);
      }
      } catch (error) {
        console.error('Camera structuralization failed:', error);
    } finally {
      setIsProcessing(false);
    }
  };

  const switchCamera = async () => {
    if (stream) {
      stream.getTracks().forEach(track => track.stop());
    }
    
    const currentFacingMode = stream?.getVideoTracks()[0]?.getSettings()?.facingMode || 'environment';
    const newFacingMode = currentFacingMode === 'environment' ? 'user' : 'environment';
    
    try {
      const mediaStream = await navigator.mediaDevices.getUserMedia({ 
        video: { facingMode: newFacingMode },
        audio: false 
      });
      
      if (videoRef.current) {
        videoRef.current.srcObject = mediaStream;
      }
      
      setStream(mediaStream);
    } catch (error) {
      console.error('Camera switch error:', error);
    }
  };

  return (
    <div className="camera-box" onClick={handleCameraClick}>
      <video 
        ref={videoRef}
        autoPlay 
        playsInline 
        muted 
        className="camera-video"
      />
      
      {!hasPermission && (
        <div className="camera-permission">{permissionMessage}</div>
      )}
      
      {/* Red outline box where user clicked */}
      {showOutline && clickPosition && (
        <div className="click-outline" style={{ left: `${clickPosition.x - 30}px`, top: `${clickPosition.y - 30}px` }}></div>
      )}
      
      {showSwitchCamera && (
        <div className="camera-switch">
          <button 
            type="button" 
            className="btn-ghost"
            onClick={(e) => {
              e.stopPropagation();
              switchCamera();
            }}
          >
            <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
              <circle cx="12" cy="12" r="5" fill="currentColor" opacity="0.8"/>
              <circle cx="12" cy="12" r="9" stroke="currentColor" fill="none"/>
            </svg>
            
          </button>
        </div>
      )}
      
      
    </div>
  );
};

const PhotoSection = ({ isProcessing, setIsProcessing, setCurrentAnalysisType, onAnalysisComplete }) => {
  const fileInputRef = useRef(null);
  const [isDragOver, setIsDragOver] = useState(false);
  const { checkPaymentRequired } = usePayment();
  const { analyzeImage } = useApi();
  const isMobile = isMobileDevice();

  const handleFileSelect = async (event) => {
    const file = event.target.files[0];
    if (!file) return;

    if (checkPaymentRequired()) {
      return;
    }

    setIsProcessing(true);
    setCurrentAnalysisType('photo');
    
    try {
      const reader = new FileReader();
      reader.onload = async (e) => {
        try {
          const result = await analyzeImage(e.target.result, file.name);
          if (onAnalysisComplete) {
            onAnalysisComplete(result);
          }
      } catch (error) {
        console.error('Photo structuralization failed:', error);
        } finally {
          setIsProcessing(false);
        }
      };
      reader.readAsDataURL(file);
    } catch (error) {
      console.error('File reading failed:', error);
      setIsProcessing(false);
    }
  };

  const handleDragOver = (e) => {
    e.preventDefault();
    setIsDragOver(true);
  };

  const handleDragLeave = (e) => {
    e.preventDefault();
    setIsDragOver(false);
  };

  const handleDrop = (e) => {
    e.preventDefault();
    setIsDragOver(false);
    
    const files = e.dataTransfer.files;
    if (files.length > 0) {
      const file = files[0];
      if (file.type.startsWith('image/')) {
        handleFileSelect({ target: { files: [file] } });
      }
    }
  };

  // Mobile: Simple button interface
  if (isMobile) {
    return (
      <div>
        <input
          ref={fileInputRef}
          type="file"
          accept="image/*"
          onChange={handleFileSelect}
          className="hidden"
        />
        <button
          onClick={() => fileInputRef.current?.click()}
          className="btn-upload"
        >
          <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" aria-hidden>
            <rect x="3" y="3" width="18" height="18" rx="2" ry="2"/>
            <circle cx="9" cy="9" r="2"/>
            <path d="M21 15l-3.086-3.086a2 2 0 0 0-2.828 0L6 21"/>
          </svg>
          
        </button>
      </div>
    );
  }

  // Desktop: Click area with red outline and drag-and-drop
  return (
    <div>
      <input
        ref={fileInputRef}
        type="file"
        accept="image/*"
        onChange={handleFileSelect}
        className="hidden"
      />
      <div
        onClick={() => fileInputRef.current?.click()}
        onDragOver={handleDragOver}
        onDragLeave={handleDragLeave}
        onDrop={handleDrop}
        className={`dropzone${isDragOver ? ' dragover' : ''}`}
      >
        <svg width="48" height="48" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.5">
          <rect x="3" y="3" width="18" height="18" rx="2" ry="2"/>
          <circle cx="9" cy="9" r="2"/>
          <path d="M21 15l-3.086-3.086a2 2 0 0 0-2.828 0L6 21"/>
        </svg>
        
      </div>
    </div>
  );
};

const LinkSection = ({ isProcessing, setIsProcessing, onAnalysisComplete }) => {
  const [imageUrl, setImageUrl] = useState('');
  const [urlError, setUrlError] = useState('');
  const { checkPaymentRequired } = usePayment();
  const { analyzeImage } = useApi();

  const validateUrl = (url) => {
    try {
      new URL(url);
      return url.match(/\.(jpeg|jpg|gif|png|webp)$/i) ? null : 'URL must point to an image file';
    } catch {
      return 'Please enter a valid URL';
    }
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    
    const validation = validateUrl(imageUrl);
    if (validation) {
      setUrlError(validation);
      return;
    }

    if (checkPaymentRequired()) {
      return;
    }

    setIsProcessing(true);
    setUrlError('');
    
    try {
      const result = await analyzeImage(imageUrl, 'Image from URL');
      if (onAnalysisComplete) {
        onAnalysisComplete(result);
      }
      setImageUrl('');
      } catch (error) {
        console.error('URL structuralization failed:', error);
        setUrlError('Failed to structuralize from URL');
    } finally {
      setIsProcessing(false);
    }
  };

  return (
    <form onSubmit={handleSubmit}>
      <div className="url-row">
        <input
          type="url"
          placeholder="Enter image URL..."
          value={imageUrl}
          onChange={(e) => {
            setImageUrl(e.target.value);
            if (urlError) setUrlError('');
          }}
          className={`input-base url-input${urlError ? ' input-error' : ''}`}
        />
        <button
          type="submit"
          disabled={!imageUrl.trim() || isProcessing}
          className="btn-icon"
        >
          →
        </button>
      </div>
      {urlError && (
        <div className="error-text">
          {urlError}
        </div>
      )}
    </form>
  );
};

// Results Display Components (Bottom Content)
 const MoleculeViewer = ({ molecularData }) => {
  const ref = useRef(null);
  const [status, setStatus] = useState('loading');

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

        // Keep existing viewer visible until new SDF is ready to avoid flash
        const hadViewer = ref.current.childNodes && ref.current.childNodes.length > 0;
        let sdfContent = null;
        if (molecularData.sdfData && molecularData.sdfData.startsWith('file://')) {
          const rawPath = molecularData.sdfData.replace('file://', '');
          const url = rawPath.startsWith('/sdf_files/')
            ? rawPath
            : `/sdf_files/${rawPath.split('/').pop()}`;
          const response = await fetch(url);
          if (response.ok) {
            sdfContent = await response.text();
            // SDF fetched
          } else {
            console.log(`SDF fetch failed for ${molecularData.name}: HTTP ${response.status}`);
          }
        }

        if (sdfContent) {
          // Now replace the canvas with a fresh viewer
          if (!ref.current) return;
          ref.current.innerHTML = '';
          const viewer = window.$3Dmol.createViewer(ref.current, {
            backgroundColor: '#000000',
            antialias: true,
            defaultcolors: window.$3Dmol.rasmolElementColors
          });
          if (typeof viewer.setBackgroundColor === 'function') {
            viewer.setBackgroundColor(0x000000, 0);
          }
          if (typeof viewer.removeAllModels === 'function') viewer.removeAllModels();
          if (typeof viewer.clear === 'function') viewer.clear();
          viewer.addModel(sdfContent, 'sdf');
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
          viewer.zoomTo();
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
        console.log(`Molecule load failed for ${molecularData.name}: ${error.message || 'Unknown error'}`);
        setStatus('failed');
      }
    };

    initialize();
    return () => {
      cancelled = true;
      if (ref.current) {
        try {
          ref.current.innerHTML = '';
        } catch (e) {
          // Ignore cleanup errors
        }
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
        {status === 'loading' && ' ⏳'} {status === 'failed' && ' ❌'}
      </div>
      <div 
        ref={ref}
        className="viewer"
      >
        {status === 'loading' && '⏳'}
        {status === 'failed' && '❌'}
      </div>
    </div>
  );
};

const MolecularColumn = ({ column, onRemove, showRemove = true }) => {
  // Render column header and viewers
  
  return (
    <div style={styles.column}>
      {/* GUARANTEED VISIBLE HEADER */}
      <div className="column-header">
        <div className="column-meta">
          <div aria-hidden role="img" className="column-icon">🧪</div>
          <div className="column-title">
            {column.query} {column.loading && '⏳'}
          </div>
        </div>
        {showRemove && (
          <button 
            onClick={onRemove}
            className="btn-ghost"
          >
            ×
          </button>
        )}
      </div>

      {/* Failure indicator */}
      {column.failed && (
        <div className="alert-warning">
          💡 AI failed to analyze - please report this
        </div>
      )}

      {column.loading && column.viewers.length === 0 && (
        <div className="analyzing">
          Analyzing molecules...
        </div>
      )}

      {column.viewers.map((mol, idx) => (
        <MoleculeViewer key={idx} molecularData={mol} />
      ))}
    </div>
  );
};

function App() {
  const [objectInput, setObjectInput] = useState('');
  const [cameraMode, setCameraMode] = useState(false);
  const [photoMode, setPhotoMode] = useState(false);
  const [linkMode, setLinkMode] = useState(false);
  const [isProcessing, setIsProcessing] = useState(false);
  const [columns, setColumns] = useState([]);
  const [error, setError] = useState('');
  const [autoVisualMode, setAutoVisualMode] = useState(process.env.NODE_ENV === 'development');
  const [showSettings, setShowSettings] = useState(false);
  const [columnMode, setColumnMode] = useState('accumulate'); // 'replace' or 'accumulate'

  const { structuresFromText: analyzeText, generateSDFs, twoNamesToSdf } = useApi();

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

  // Handle keyboard shortcuts
  useEffect(() => {
    const handleKeyDown = (event) => {
      if (event.target.tagName === 'INPUT' || event.target.tagName === 'TEXTAREA') {
        return;
      }

      const isMac = navigator.userAgent.toUpperCase().indexOf('MAC') >= 0;
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

    
    // Create column based on mode setting
    const columnId = Date.now();
    const newColumn = {
      id: columnId,
      query: value,
      viewers: [],
      loading: true,
      failed: false
    };
    
    if (columnMode === 'replace') {
      setColumns([newColumn]); // Replace existing
    } else {
      setColumns(prev => [...prev, newColumn]); // Add to existing
    }
    
    try {
      // Use two-names→SDF flow only when EXACTLY two names are provided
      const parts = value.split(',').map(s => s.trim()).filter(Boolean);
      let result;
      let usedTwoNames = false;
      if (parts.length === 2) {
        try {
          const resp = await twoNamesToSdf(parts, false);
          result = { object: value, molecules: resp.molecules || [] };
          usedTwoNames = true;
        } catch (e) {
          // Fallback to standard analyzer if two-names call errors
          result = await analyzeText(value);
        }
      } else {
        result = await analyzeText(value);
      }
      const molecules = result.molecules || result.chemicals || [];

      if (molecules && molecules.length > 0) {
        // If two-names was used but we received no SDFs and no SMILES, fallback to analyzer
        if (usedTwoNames) {
          const hasAnySdf = molecules.some(m => !!m.sdfPath);
          const hasAnySmiles = molecules.some(m => !!m.smiles);
          if (!hasAnySdf && !hasAnySmiles) {
            try {
              const fallback = await analyzeText(value);
              result = fallback;
            } catch (_) {
              // keep original result
            }
          }
        }
        // Refresh molecules after possible fallback
        const molecules = result.molecules || result.chemicals || [];
        // Prefer precomputed SDFs with canonical names
        const precomputed = molecules.filter(m => m.sdfPath);
        if (precomputed.length > 0) {
          const viewers = precomputed.map(m => ({
            name: m.name || value,
            sdfData: m.sdfPath.startsWith('file://') ? m.sdfPath : `file://${m.sdfPath}`,
            smiles: m.smiles
          }));
          setColumns(prev => prev.map(col => (
            col.id === columnId ? { ...col, viewers, loading: false, failed: false } : col
          )));
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
              setColumns(prev => prev.map(col => (
                col.id === columnId ? { ...col, viewers, loading: false, failed: false } : col
              )));
            } catch (sdfError) {
              console.error('SDF generation failed:', sdfError);
              setError('Failed to generate molecular structures');
              setColumns(prev => prev.map(col => (
                col.id === columnId ? { ...col, loading: false, failed: true } : col
              )));
            }
          } else {
            setError('No valid molecules returned by the analyzer.');
            setColumns(prev => prev.map(col => (
              col.id === columnId ? { ...col, loading: false, failed: true } : col
            )));
          }
        }
      } else {
        setError('No molecules found for this input.');
        setColumns(prev => prev.map(col => (
          col.id === columnId ? { ...col, loading: false, failed: true } : col
        )));
      }
      
      setObjectInput('');
    } catch (error) {
      console.error('Analysis failed:', error);
      setError('Analysis failed. Please try again.');
      // Mark column as failed
      setColumns(prev => prev.map(col => (
        col.id === columnId ? { ...col, loading: false, failed: true } : col
      )));
    } finally {
      setIsProcessing(false);
    }
  }, [isProcessing, analyzeText, generateSDFs, columnMode]);

  const handleAnalysisComplete = useCallback(async (result) => {
    const molecules = result?.molecules || result?.chemicals || [];
    const objectLabel = (result && (result.object || (result.result && result.result.object))) || (cameraMode ? 'Camera capture' : 'Image capture');

    // Create column based on mode setting
    const columnId = Date.now();
    const newColumn = {
      id: columnId,
      query: objectLabel,
      viewers: [],
      loading: true,
      failed: false
    };
    
    if (columnMode === 'replace') {
      setColumns([newColumn]); // Replace existing
    } else {
      setColumns(prev => [...prev, newColumn]); // Add to existing
    }

    if (molecules && molecules.length > 0) {
      const smilesArray = molecules.map(mol => mol.smiles).filter(Boolean);

      if (smilesArray.length > 0) {
        try {
          const sdfResult = await generateSDFs(smilesArray, false);

          const smilesToSdf = new Map();
          sdfResult.sdfPaths?.forEach((p, i) => { smilesToSdf.set(smilesArray[i], p); });

          const viewers = molecules.map((mol) => {
            const sdfPath = smilesToSdf.get(mol.smiles);
            return {
              name: mol.name || (mol.smiles && SMILES_NAME_MAP[mol.smiles]) || mol.smiles || objectLabel,
              sdfData: sdfPath ? `file://${sdfPath}` : null,
              smiles: mol.smiles
            };
          });

          // Update the column with viewers
          setColumns(prev => prev.map(col => (
            col.id === columnId ? { ...col, viewers, loading: false, failed: false } : col
          )));
        } catch (sdfError) {
          console.error('SDF generation failed:', sdfError);
          setColumns(prev => prev.map(col => (
            col.id === columnId ? { ...col, loading: false, failed: true } : col
          )));
        }
      } else {
        setColumns(prev => prev.map(col => (
          col.id === columnId ? { ...col, loading: false, failed: true } : col
        )));
      }
    } else {
      setColumns(prev => prev.map(col => (
        col.id === columnId ? { ...col, loading: false, failed: true } : col
      )));
    }

    setError('');
  }, [generateSDFs, cameraMode, columnMode]);

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
    <PaymentProvider config={PAYMENT_CONFIG}>
      <div className="app">
        {/* Settings gear icon in top right - OUTSIDE main to avoid conflicts */}
        <button
          onClick={() => {
            console.log('Settings button clicked, current state:', showSettings);
            setShowSettings(!showSettings);
          }}
          className="settings-btn"
          title="Settings"
        >
          ⚙️
        </button>
        
        {/* Settings Modal - OUTSIDE main for proper positioning */}
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
                  console.log('Column mode changed to:', e.target.value);
                  setColumnMode(e.target.value);
                }}
                className="select"
              >
                <option value="replace" style={{ background: '#222', color: '#fff' }}>
                  Replace (Default) - New analysis displaces previous
                </option>
                <option value="accumulate" style={{ background: '#222', color: '#fff' }}>
                  Accumulate - Add columns side by side
                </option>
              </select>
            </div>

            <button
              onClick={() => setShowSettings(false)}
              className="settings-close"
            >
              Close
            </button>
          </div>
        )}
        
        <div className="main">          
          <div className="input-section">
            <div className="mb-15">
              <TextInput 
                value={objectInput}
                onChange={setObjectInput}
                onSubmit={handleTextAnalysis}
                isProcessing={isProcessing}
                error={error}
              />
            </div>

            <div className="mb-15">
              <ModeSelector
                cameraMode={cameraMode}
                setCameraMode={setCameraMode}
                photoMode={photoMode}
                setPhotoMode={setPhotoMode}
                linkMode={linkMode}
                setLinkMode={setLinkMode}
              />
            </div>

            {cameraMode && (
              <div className="mb-15">
                <CameraSection
                  isProcessing={isProcessing}
                  setIsProcessing={setIsProcessing}
                  setCurrentAnalysisType={() => {}}
                  onAnalysisComplete={handleAnalysisComplete}
                />
              </div>
            )}

            {photoMode && (
              <div className="mb-15">
                <PhotoSection
                  isProcessing={isProcessing}
                  setIsProcessing={setIsProcessing}
                  setCurrentAnalysisType={() => {}}
                  onAnalysisComplete={handleAnalysisComplete}
                />
              </div>
            )}

            {linkMode && (
              <div className="mb-15">
                <LinkSection
                  isProcessing={isProcessing}
                  setIsProcessing={setIsProcessing}
                  onAnalysisComplete={handleAnalysisComplete}
                />
              </div>
            )}
          </div>



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
    </PaymentProvider>
  );
}

export default App;