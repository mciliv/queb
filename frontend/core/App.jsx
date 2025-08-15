import React, { useState, useEffect, useCallback, useRef } from 'react';
import { PaymentProvider, usePayment } from '../components/ui/PaymentContext';
import { useApi } from '../hooks/useApi';
import { PRESET_VISUAL_TESTS, TEST_MOLECULES, SMILES_NAME_MAP } from './constants.js';
import '../assets/style.css';

// Device detection utility
const isMobileDevice = () => {
  return /Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini/i.test(navigator.userAgent) ||
         (window.innerWidth <= 768 && 'ontouchstart' in window);
};

// Configuration
const PAYMENT_CONFIG = {
  enabled: false, // Set to true to enable payment functionality
  devMode: window.location.hostname === 'localhost' || window.location.hostname === '127.0.0.1'
};

// Styles
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
    letterSpacing: '0.01em',
    userSelect: 'text'
  },
  mainLayout: {
    padding: '20px',
    userSelect: 'text'
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
    marginBottom: '30px',
    maxWidth: '600px'
  },
  columnsContainer: {
    display: 'flex',
    flexDirection: 'row',
    gap: '20px',
    overflowX: 'auto',
    position: 'relative',
    zIndex: 1,
    userSelect: 'text'
  },
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
    <div style={{ position: 'relative', width: '100%', marginBottom: '20px' }}>
      <div style={{ position: 'relative', display: 'flex', alignItems: 'center', gap: '8px' }}>
        <input
          ref={inputRef}
          id="object-input"
          type="text"
          placeholder="Type to analyze molecules..."
          style={{
            width: '100%',
            padding: '12px 16px',
            background: 'rgba(255, 255, 255, 0.08)',
            border: displayError ? '1px solid rgba(255, 100, 100, 0.5)' : 'none',
            borderRadius: '8px',
            color: '#ffffff',
            fontSize: '14px',
            outline: 'none',
            WebkitAppearance: 'none',
            appearance: 'none',
            fontFamily: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif'
          }}
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
          <div style={{
            position: 'absolute',
            right: '12px',
            top: '50%',
            transform: 'translateY(-50%)',
            background: 'rgba(255, 255, 255, 0.1)',
            border: '1px solid rgba(255, 255, 255, 0.2)',
            borderRadius: '4px',
            padding: '4px 8px',
            fontSize: '11px',
            color: 'rgba(255, 255, 255, 0.6)',
            fontFamily: 'monospace',
            pointerEvents: 'none',
            userSelect: 'none'
          }}>
            {keyboardHint}
          </div>
        )}
      </div>
      {displayError && (
        <div id="input-error" style={{ color: '#ff6b6b', fontSize: '12px', marginTop: '8px', padding: '0 4px' }} role="alert">
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
    <div style={{ display: 'flex', gap: '10px', marginTop: '10px' }}>
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
    console.log('Camera clicked at:', { x, y });
    setClickPosition({ x, y });
    setShowOutline(true);
    
    // Hide outline after 2 seconds
    setTimeout(() => {
      console.log('Hiding red outline box');
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
    <div style={{
      position: 'relative',
      width: '100%',
      maxWidth: '400px',
      borderRadius: '8px',
      overflow: 'hidden',
      cursor: 'crosshair',
      background: '#000',
      aspectRatio: '4/3'
    }} onClick={handleCameraClick}>
      <video 
        ref={videoRef}
        autoPlay 
        playsInline 
        muted 
        style={{ width: '100%', height: '100%', objectFit: 'cover' }}
      />
      
      {!hasPermission && (
        <div style={{
          position: 'absolute',
          top: '50%',
          left: '50%',
          transform: 'translate(-50%, -50%)',
          color: '#ffffff',
          textAlign: 'center',
          padding: '20px',
          background: 'rgba(0, 0, 0, 0.8)',
          borderRadius: '8px',
          fontSize: '14px'
        }}>{permissionMessage}</div>
      )}
      
      {/* Red outline box where user clicked */}
      {showOutline && clickPosition && (
        <div style={{
          position: 'absolute',
          left: `${clickPosition.x - 30}px`,
          top: `${clickPosition.y - 30}px`,
          width: '60px',
          height: '60px',
          border: '3px solid #ff3333',
          borderRadius: '8px',
          pointerEvents: 'none',
          animation: 'pulse 0.5s ease-in-out',
          zIndex: 100,
          boxShadow: '0 0 10px rgba(255, 51, 51, 0.5)',
          background: 'rgba(255, 51, 51, 0.1)'
        }}></div>
      )}
      
      {showSwitchCamera && (
        <div style={{ position: 'absolute', top: '16px', right: '16px' }}>
          <button 
            type="button" 
            style={{
              display: 'flex',
              alignItems: 'center',
              gap: '8px',
              background: 'rgba(0, 0, 0, 0.7)',
              color: '#ffffff',
              border: 'none',
              borderRadius: '20px',
              padding: '8px 12px',
              fontSize: '12px',
              cursor: 'pointer'
            }}
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
          style={{ display: 'none' }}
        />
        <button
          onClick={() => fileInputRef.current?.click()}
          style={{
            width: '100%',
            padding: '12px',
            background: 'rgba(255, 255, 255, 0.08)',
            border: 'none',
            borderRadius: '8px',
            color: '#ffffff',
            fontSize: '14px',
            cursor: 'pointer',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            gap: '8px'
          }}
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
        style={{ display: 'none' }}
      />
      <div
        onClick={() => fileInputRef.current?.click()}
        onDragOver={handleDragOver}
        onDragLeave={handleDragLeave}
        onDrop={handleDrop}
        style={{
          width: '100%',
          height: '200px',
          border: 'none',
          borderRadius: '8px',
          display: 'flex',
          flexDirection: 'column',
          alignItems: 'center',
          justifyContent: 'center',
          cursor: 'pointer',
          background: isDragOver ? 'rgba(255, 255, 255, 0.08)' : 'rgba(255, 255, 255, 0.05)',
          transition: 'all 0.2s ease',
          gap: '12px'
        }}
      >
        <svg width="48" height="48" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.5" style={{ opacity: 0.6 }}>
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
      <div style={{ display: 'flex', gap: '8px', alignItems: 'center' }}>
        <input
          type="url"
          placeholder="Enter image URL..."
          value={imageUrl}
          onChange={(e) => {
            setImageUrl(e.target.value);
            if (urlError) setUrlError('');
          }}
          style={{
            flex: 1,
            padding: '12px 16px',
            background: 'rgba(255, 255, 255, 0.08)',
            border: urlError ? '1px solid rgba(255, 100, 100, 0.5)' : 'none',
            borderRadius: '8px',
            color: '#ffffff',
            fontSize: '14px',
            outline: 'none'
          }}
        />
        <button
          type="submit"
          disabled={!imageUrl.trim() || isProcessing}
          style={{
            background: 'rgba(255, 255, 255, 0.1)',
            border: 'none',
            borderRadius: '6px',
            color: '#ffffff',
            padding: '8px 12px',
            cursor: imageUrl.trim() && !isProcessing ? 'pointer' : 'not-allowed',
            fontSize: '16px',
            minWidth: '40px',
            height: '40px',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            opacity: imageUrl.trim() && !isProcessing ? 1 : 0.5
          }}
        >
          →
        </button>
      </div>
      {urlError && (
        <div style={{ color: '#ff6b6b', fontSize: '12px', marginTop: '8px' }}>
          {urlError}
        </div>
      )}
    </form>
  );
};

// Results Display Components (Bottom Content)
const SimpleMoleculeViewer = ({ molecularData }) => {
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

        console.log(`Loading molecule: ${molecularData.name} with SMILES: ${molecularData.smiles}`);
        
        const viewer = window.$3Dmol.createViewer(ref.current, {
          backgroundColor: '#000000',
          antialias: true,
          defaultcolors: window.$3Dmol.rasmolElementColors
        });
        // Make background transparent while using a valid color (avoid 'transparent' error)
        if (typeof viewer.setBackgroundColor === 'function') {
          viewer.setBackgroundColor(0x000000, 0);
        }

        let sdfContent = null;
        if (molecularData.sdfData && molecularData.sdfData.startsWith('file://')) {
          const path = molecularData.sdfData.replace('file://', '');
          const response = await fetch(`/sdf_files/${path.split('/').pop()}`);
          if (response.ok) {
            sdfContent = await response.text();
            console.log(`SDF loaded for ${molecularData.name}`);
          } else {
            console.error(`Failed to fetch SDF for ${molecularData.name}`);
          }
        }

        if (sdfContent) {
          // Ensure no residual models/styles from prior renders
          if (typeof viewer.removeAllModels === 'function') {
            viewer.removeAllModels();
          }
          if (typeof viewer.clear === 'function') {
            viewer.clear();
          }
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
          console.log(`${molecularData.name} rendered successfully with element colors`);
        } else {
          setStatus('failed');
        }
      } catch (error) {
        console.error(`Error loading ${molecularData.name}:`, error);
        setStatus('failed');
      }
    };

    initialize();
    return () => {
      cancelled = true;
      if (ref.current) {
        ref.current.innerHTML = '';
      }
    };
  }, [molecularData.sdfData, molecularData.name, molecularData.smiles]);

  return (
    <div style={{ marginBottom: '12px' }}>
      <div style={{ 
        background: 'transparent', 
        padding: '8px', 
        marginBottom: '5px',
        borderRadius: '4px',
        fontSize: '14px',
        userSelect: 'text'
      }}>
        {molecularData.name && (
          <a 
            href={`https://en.wikipedia.org/wiki/${encodeURIComponent(molecularData.name)}`}
            target="_blank"
            rel="noopener noreferrer"
            style={{
              color: '#ffffff',
              textDecoration: 'none'
            }}
          >
            {molecularData.name}
          </a>
        )}
        {status === 'loading' && ' ⏳'} {status === 'failed' && ' ❌'}
      </div>
      <div 
        ref={ref}
        style={{ 
          height: '200px', 
          width: '100%', 
          background: 'transparent', 
          borderRadius: '4px',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          color: 'rgba(255,255,255,0.5)',
          position: 'relative',
          overflow: 'hidden',
          zIndex: 1
        }}
      >
        {status === 'loading' && '⏳'}
        {status === 'failed' && '❌'}
      </div>
    </div>
  );
};

const MolecularColumn = ({ column, onRemove, showRemove = true }) => {
  // Debug: Log column data to ensure headers have content
  console.log('Rendering column header:', { query: column.query, loading: column.loading, failed: column.failed });
  
  return (
    <div style={styles.column}>
      {/* GUARANTEED VISIBLE HEADER */}
      <div style={{ 
        display: 'flex', 
        justifyContent: 'space-between', 
        alignItems: 'center', 
        marginBottom: '20px', 
        userSelect: 'text',
        minHeight: '24px',
        position: 'relative',
        zIndex: 10,
        visibility: 'visible',
        opacity: 1,
        background: 'rgba(255, 255, 255, 0.02)', // Subtle background to ensure visibility
        padding: '8px 0',
        borderRadius: '4px'
      }}>
        <div style={{ 
          display: 'flex', 
          alignItems: 'center', 
          gap: '8px', 
          userSelect: 'text',
          flex: 1
        }}>
          <div aria-hidden role="img" style={{ 
            opacity: 0.8, 
            fontSize: '16px', 
            userSelect: 'none',
            minWidth: '16px'
          }}>🧪</div>
          <div style={{ 
            fontSize: '14px', 
            opacity: 0.7, 
            userSelect: 'text',
            color: '#ffffff',
            fontWeight: '400'
          }}>
            {column.query} {column.loading && '⏳'}
          </div>
        </div>
        {showRemove && (
          <button 
            onClick={onRemove}
            style={{
              background: 'transparent',
              border: 'none',
              color: '#ffffff',
              fontSize: '20px',
              cursor: 'pointer',
              userSelect: 'none',
              minWidth: '24px',
              minHeight: '24px',
              padding: '2px',
              flexShrink: 0
            }}
          >
            ×
          </button>
        )}
      </div>

      {/* Failure indicator */}
      {column.failed && (
        <div style={{
          marginBottom: '16px',
          padding: '10px 12px',
          background: 'rgba(255, 193, 7, 0.1)',
          border: '1px solid rgba(255, 193, 7, 0.3)',
          borderRadius: '4px',
          fontSize: '13px',
          color: 'rgba(255, 193, 7, 0.9)'
        }}>
          💡 AI failed to analyze - please report this
        </div>
      )}

      {column.loading && column.viewers.length === 0 && (
        <div style={{ 
          padding: '20px', 
          textAlign: 'center', 
          color: 'rgba(255,255,255,0.5)',
          fontSize: '13px',
          userSelect: 'text'
        }}>
          Analyzing molecules...
        </div>
      )}

      {column.viewers.map((mol, idx) => (
        <SimpleMoleculeViewer key={idx} molecularData={mol} />
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

  const { analyzeText, generateSDFs } = useApi();

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
    if (PAYMENT_CONFIG.devMode) {
      console.log('🔧 Auto-enabling developer mode for localhost');
    }
    console.log('✅ Molecular analysis app initialized');
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
      const result = await analyzeText(value);
      const molecules = result.molecules || result.chemicals || [];
      
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
                name: mol.name || (mol.smiles && SMILES_NAME_MAP[mol.smiles]) || mol.smiles || value,
                sdfData: sdfPath ? `file://${sdfPath}` : null,
                smiles: mol.smiles
              };
            });
            
            // Update the newly created column by id
            setColumns(prev => prev.map(col => (
              col.id === columnId ? { ...col, viewers, loading: false, failed: false } : col
            )));
          } catch (sdfError) {
            console.error('SDF generation failed:', sdfError);
            setError('Failed to generate molecular structures');
            // Mark this column as failed
            setColumns(prev => prev.map(col => (
              col.id === columnId ? { ...col, loading: false, failed: true } : col
            )));
          }
        } else {
          setError('No valid SMILES found in the analysis results.');
          setColumns(prev => prev.map(col => (
            col.id === columnId ? { ...col, loading: false, failed: true } : col
          )));
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
      await new Promise(resolve => setTimeout(resolve, 500));
      
      for (const test of PRESET_VISUAL_TESTS) {
        if (cancelled) break;
        try {
          // Generate SDFs for multiple compounds per object
          const smilesArray = test.smilesList || [];
          if (smilesArray.length === 0) continue;
          const sdfResult = await generateSDFs(smilesArray, false);
          const viewers = smilesArray.map((sm, idx) => ({
            name: SMILES_NAME_MAP[sm] || sm,
            sdfData: sdfResult.sdfPaths?.[idx] ? `file://${sdfResult.sdfPaths[idx]}` : null,
            smiles: sm
          }));
          setColumns(prev => ([...prev, { id: Date.now() + Math.random(), query: test.label, viewers }]));
        } catch (e) {
          console.log(`Failed to create column for ${test.label}:`, e);
        }
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
          style={{
            position: 'fixed',
            top: '20px',
            right: '20px',
            width: '44px',
            height: '44px',
            borderRadius: '22px',
            background: 'rgba(255, 255, 255, 0.15)',
            border: '1px solid rgba(255, 255, 255, 0.2)',
            color: '#ffffff',
            cursor: 'pointer',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            fontSize: '18px',
            zIndex: 9999,
            opacity: 1,
            visibility: 'visible',
            boxShadow: '0 4px 12px rgba(0, 0, 0, 0.3)'
          }}
          title="Settings"
        >
          ⚙️
        </button>
        
        {/* Settings Modal - OUTSIDE main for proper positioning */}
        {showSettings && (
          <div style={{
            position: 'fixed',
            top: '75px',
            right: '20px',
            background: 'rgba(20, 20, 20, 0.95)',
            border: '1px solid rgba(255, 255, 255, 0.2)',
            borderRadius: '8px',
            padding: '20px',
            minWidth: '280px',
            zIndex: 9998,
            backdropFilter: 'blur(10px)',
            boxShadow: '0 8px 24px rgba(0, 0, 0, 0.4)'
          }}>
            <h3 style={{ margin: '0 0 15px 0', fontSize: '16px', color: '#ffffff' }}>Settings</h3>
            
            <div style={{ marginBottom: '15px' }}>
              <label style={{ fontSize: '14px', color: '#ffffff', marginBottom: '8px', display: 'block' }}>
                Column Behavior:
              </label>
              <select
                value={columnMode}
                onChange={(e) => {
                  console.log('Column mode changed to:', e.target.value);
                  setColumnMode(e.target.value);
                }}
                style={{
                  width: '100%',
                  padding: '8px 12px',
                  background: 'rgba(255, 255, 255, 0.1)',
                  border: '1px solid rgba(255, 255, 255, 0.2)',
                  borderRadius: '4px',
                  color: '#ffffff',
                  fontSize: '14px'
                }}
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
              style={{
                background: 'rgba(255, 255, 255, 0.1)',
                border: '1px solid rgba(255, 255, 255, 0.2)',
                borderRadius: '4px',
                color: '#ffffff',
                padding: '8px 16px',
                cursor: 'pointer',
                fontSize: '14px',
                width: '100%'
              }}
            >
              Close
            </button>
          </div>
        )}
        
        <div className="main">          
          <div style={styles.inputSection}>
            <div style={{ marginBottom: '15px' }}>
              <TextInput 
                value={objectInput}
                onChange={setObjectInput}
                onSubmit={handleTextAnalysis}
                isProcessing={isProcessing}
                error={error}
              />
            </div>

            <div style={{ marginBottom: '15px' }}>
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
              <div style={{ marginBottom: '15px' }}>
                <CameraSection
                  isProcessing={isProcessing}
                  setIsProcessing={setIsProcessing}
                  setCurrentAnalysisType={() => {}}
                  onAnalysisComplete={handleAnalysisComplete}
                />
              </div>
            )}

            {photoMode && (
              <div style={{ marginBottom: '15px' }}>
                <PhotoSection
                  isProcessing={isProcessing}
                  setIsProcessing={setIsProcessing}
                  setCurrentAnalysisType={() => {}}
                  onAnalysisComplete={handleAnalysisComplete}
                />
              </div>
            )}

            {linkMode && (
              <div style={{ marginBottom: '15px' }}>
                <LinkSection
                  isProcessing={isProcessing}
                  setIsProcessing={setIsProcessing}
                  onAnalysisComplete={handleAnalysisComplete}
                />
              </div>
            )}
          </div>



          <div className="columns" style={{ display: 'flex', flexDirection: 'row', gap: '20px', overflowX: 'auto' }}>
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