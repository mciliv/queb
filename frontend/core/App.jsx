import React, { useState, useEffect, useCallback, useRef } from 'react';
import { PaymentProvider, usePayment } from '../components/ui/PaymentContext';
import { useApi } from '../hooks/useApi';
import { PRESET_VISUAL_TESTS, TEST_MOLECULES, SMILES_NAME_MAP } from './constants.js';
import { PAYMENT_CONFIG, VALIDATION_PATTERNS } from '../utils/config-loader.js';
import logger from './logger.js';
import { createKeyboardHandler } from './keyboard-shortcuts.js';
import '../assets/style.css';

const isMobileDevice = () => {
  return /Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini/i.test(navigator.userAgent) ||
         (window.innerWidth <= 768 && 'ontouchstart' in window);
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
    
    const invalidPatterns = Object.values(VALIDATION_PATTERNS);
    
    if (invalidPatterns.some(pattern => pattern.test(trimmed))) {
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
          placeholder="Specify object..."
          className={`input-base${displayError ? ' input-error' : ''}`}
          value={value}
          onChange={(e) => onChange(e.target.value)}
          onKeyDown={handleKeyDown}
          aria-describedby={displayError ? 'input-error' : undefined}
        />
        {!isMobileDevice() && !value.trim() && (
          <div className="kbd-hint kbd-inside">{keyboardHint}</div>
        )}
        
        {value.trim() && (
          <button 
            id="object-submit"
            className="btn-icon"
            onClick={handleSubmit}
            aria-label="Structuralize"
          >
            →
         </button>
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
  const isMobile = isMobileDevice();
  
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



  return (
    <div className="mode-row">
      <button 
        className={`mode-btn${cameraMode ? ' active' : ''}`}
        onClick={() => handleModeSelect('camera')}
        title={isMobile ? "Capture from camera" : "Capture from camera (⌥C)"}
      >
        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
          <circle cx="12" cy="12" r="5" fill="currentColor" opacity="0.8"/>
          <circle cx="12" cy="12" r="9" stroke="currentColor" fill="none"/>
        </svg>
        {!isMobile && (
          <span className="mode-btn-shortcut">
            {navigator.userAgent.toUpperCase().indexOf('MAC') >= 0 ? '⌘M' : 'Ctrl+M'}
          </span>
        )}
      </button>

      <button 
        className={`mode-btn${photoMode ? ' active' : ''}`}
        onClick={() => handleModeSelect('photo')}
        title={isMobile ? "Upload image" : "Upload image (⌥P)"}
      >
        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
          <rect x="3" y="3" width="18" height="18" rx="2" ry="2"/>
          <circle cx="9" cy="9" r="2"/>
          <path d="M21 15l-3.086-3.086a2 2 0 0 0-2.828 0L6 21"/>
        </svg>
        {!isMobile && (
          <span className="mode-btn-shortcut">
            {navigator.userAgent.toUpperCase().indexOf('MAC') >= 0 ? '⌘P' : 'Ctrl+P'}
          </span>
        )}
      </button>

      <button 
        className={`mode-btn${linkMode ? ' active' : ''}`}
        onClick={() => handleModeSelect('link')}
        title={isMobile ? "Enter image link" : "Enter image link (⌥L)"}
      >
        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
          <path d="M10 13a5 5 0 0 0 7.54.54l3-3a5 5 0 0 0-7.07-7.07l-1.72 1.71"/>
          <path d="M14 11a5 5 0 0 0-7.54-.54l-3 3a5 5 0 0 0 7.07 7.07l1.71-1.71"/>
        </svg>
        {!isMobile && (
          <span className="mode-btn-shortcut">
            {navigator.userAgent.toUpperCase().indexOf('MAC') >= 0 ? '⌘L' : 'Ctrl+L'}
          </span>
        )}
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
  const outlineTimeoutRef = useRef(null);
  const { checkPaymentRequired } = usePayment();
  const { analyzeImage } = useApi();

  useEffect(() => {
    requestCameraAccess();
    return () => {
      if (stream) {
        stream.getTracks().forEach(track => track.stop());
      }
      if (outlineTimeoutRef.current) {
        try { clearTimeout(outlineTimeoutRef.current); } catch (_) {}
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
      logger.error('Camera access error:', error);
      setPermissionMessage('Camera access required to structuralize from camera');
      setHasPermission(false);
    }
  };

  const handleCameraClick = async (e) => {
    if (!hasPermission) {
      await requestCameraAccess();
      return;
    }
    if (!videoRef.current || isProcessing) return;

    if (checkPaymentRequired()) {
      return;
    }

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

      // Determine click position relative to video frame (accounting for object-fit: cover)
      const container = e.currentTarget;
      const rect = container.getBoundingClientRect();
      const clickX = e.clientX - rect.left;
      const clickY = e.clientY - rect.top;

      const videoAspect = width / height;
      const containerAspect = rect.width / rect.height;
      let scale = 1;
      let offsetX = 0;
      let offsetY = 0;
      if (containerAspect > videoAspect) {
        // Container is wider; height overflows
        scale = rect.width / width;
        const scaledHeight = rect.width / videoAspect;
        offsetY = (rect.height - scaledHeight) / 2;
      } else {
        // Container is taller; width overflows
        scale = rect.height / height;
        const scaledWidth = rect.height * videoAspect;
        offsetX = (rect.width - scaledWidth) / 2;
      }

      // Map click to video coordinates
      const videoClickX = (clickX - offsetX) / scale;
      const videoClickY = (clickY - offsetY) / scale;

      // Crop square side (25% of min dimension)
      const cropSide = Math.floor(Math.min(width, height) * 0.25);
      const half = Math.floor(cropSide / 2);
      // Clamp center so crop stays within frame
      const centerX = Math.min(Math.max(videoClickX, half), width - half);
      const centerY = Math.min(Math.max(videoClickY, half), height - half);
      const cropX = Math.max(0, Math.floor(centerX - half));
      const cropY = Math.max(0, Math.floor(centerY - half));

      ctx.drawImage(video, 0, 0);
      const fullImageDataUrl = canvas.toDataURL('image/jpeg', 0.8);

      // Extract cropped region as separate base64 to boost accuracy
      const cropCanvas = document.createElement('canvas');
      cropCanvas.width = cropSide;
      cropCanvas.height = cropSide;
      const cropCtx = cropCanvas.getContext('2d');
      cropCtx.drawImage(canvas, cropX, cropY, cropSide, cropSide, 0, 0, cropSide, cropSide);
      const croppedImageDataUrl = cropCanvas.toDataURL('image/jpeg', 0.8);
      
      // Show transient outline at the clicked crop region (in container coords)
      try {
        const outlineSize = cropSide * scale;
        const displayCenterX = offsetX + centerX * scale;
        const displayCenterY = offsetY + centerY * scale;
        setClickPosition({
          x: Math.round(displayCenterX - outlineSize / 2),
          y: Math.round(displayCenterY - outlineSize / 2),
          size: Math.round(outlineSize)
        });
        setShowOutline(true);
        if (outlineTimeoutRef.current) {
          try { clearTimeout(outlineTimeoutRef.current); } catch (_) {}
        }
        outlineTimeoutRef.current = setTimeout(() => setShowOutline(false), 1000);
      } catch (_) {}

      // AB test: 'coords' uses full image + coordinates; 'crop' uses cropped region
      const abVariant = Math.random() < 0.5 ? 'coords' : 'crop';
      let result;
      if (abVariant === 'coords') {
        logger.debug('Camera click payload summary', {
          variant: 'coords',
          frame: { width, height },
          mappedClick: { x: Math.round(centerX), y: Math.round(centerY) },
          crop: { cx: cropX + Math.floor(cropSide / 2), cy: cropY + Math.floor(cropSide / 2), size: cropSide },
          display: { rectWidth: Math.round(rect.width), rectHeight: Math.round(rect.height), offsetX: Math.round(offsetX), offsetY: Math.round(offsetY), scale: Number(scale.toFixed(3)) }
        });
        result = await analyzeImage(
          fullImageDataUrl,
          Math.round(centerX),
          Math.round(centerY),
          cropX + Math.floor(cropSide / 2),
          cropY + Math.floor(cropSide / 2),
          cropSide
        );
      } else {
        // For cropped case, provide center relative to the crop
        const mid = Math.floor(cropSide / 2);
        logger.debug('Camera click payload summary', {
          variant: 'crop',
          frame: { width, height },
          cropCenter: { mid },
          crop: { size: cropSide }
        });
        result = await analyzeImage(
          croppedImageDataUrl,
          mid,
          mid,
          mid,
          mid,
          cropSide
        );
      }
      if (onAnalysisComplete) {
        onAnalysisComplete({ ...result, abVariant });
      }
      } catch (error) {
        logger.error('Camera structuralization failed:', error);
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
      logger.error('Camera switch error:', error);
    }
  };

  return (
    <div
      className="camera-box"
      onClick={handleCameraClick}
    >
      {/* Visible camera preview */}
      <video
        ref={videoRef}
        autoPlay
        playsInline
        muted
        className="camera-video"
      />
      {showOutline && clickPosition && (
        <div
          className="click-outline"
          style={{
            left: clickPosition.x,
            top: clickPosition.y,
            width: clickPosition.size,
            height: clickPosition.size
          }}
        />
      )}
    </div>
  );
};

const PhotoSection = ({ isProcessing, setIsProcessing, setCurrentAnalysisType, onAnalysisComplete }) => {
  const fileInputRef = useRef(null);
  const [isDragOver, setIsDragOver] = useState(false);
  const [previewUrl, setPreviewUrl] = useState('');
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
    // Show immediate preview of selected image
    try {
      const url = URL.createObjectURL(file);
      setPreviewUrl(url);
    } catch (_) {}
    
    try {
      const reader = new FileReader();
      reader.onload = async (e) => {
        try {
          const dataUrl = e.target.result;
          // Build cropped header preview (center square, 25% of min dim)
          const img = new Image();
          img.onload = async () => {
            try {
              const w = img.naturalWidth || img.width;
              const h = img.naturalHeight || img.height;
              const cropSide = Math.floor(Math.min(w, h) * 0.25);
              const cropX = Math.max(0, Math.floor((w - cropSide) / 2));
              const cropY = Math.max(0, Math.floor((h - cropSide) / 2));
              const canvas = document.createElement('canvas');
              canvas.width = cropSide;
              canvas.height = cropSide;
              const ctx = canvas.getContext('2d');
              ctx.drawImage(img, cropX, cropY, cropSide, cropSide, 0, 0, cropSide, cropSide);
              const croppedHeaderImage = canvas.toDataURL('image/jpeg', 0.8);

              logger.debug('Photo analyze payload summary', {
                variant: 'file',
                image: { width: w, height: h },
                crop: { size: cropSide, cx: cropX + Math.floor(cropSide/2), cy: cropY + Math.floor(cropSide/2) }
              });
              const result = await analyzeImage(dataUrl);
              if (onAnalysisComplete) {
                onAnalysisComplete(result);
              }
            } catch (innerErr) {
              logger.error('Photo header crop failed:', innerErr);
            } finally {
              setIsProcessing(false);
            }
          };
          img.onerror = () => {
            // Fallback: pass original if crop fails to load
            logger.debug('Photo analyze payload summary', { variant: 'file-fallback' });
            analyzeImage(dataUrl)
              .then((result) => {
                if (onAnalysisComplete) {
                  onAnalysisComplete(result);
                }
              })
              .catch((err) => logger.error('Photo structuralization failed:', err))
              .finally(() => setIsProcessing(false));
          };
          img.src = dataUrl;
        } catch (error) {
          logger.error('File reading or processing failed:', error);
          setIsProcessing(false);
        }
      };
      reader.readAsDataURL(file);
    } catch (error) {
      logger.error('File reading failed:', error);
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
        {previewUrl && (
          <img
            src={previewUrl}
            alt="Selected"
            className="preview-image"
            onLoad={() => {
              // Revoke to release memory after first paint
              try { URL.revokeObjectURL(previewUrl); } catch (_) {}
            }}
          />
        )}
      </div>
    </div>
  );
};

const LinkSection = ({ isProcessing, setIsProcessing, onAnalysisComplete }) => {
  const [imageUrl, setImageUrl] = useState('');
  const [urlError, setUrlError] = useState('');
  const [isDragOver, setIsDragOver] = useState(false);
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
      logger.debug('Link analyze payload summary', { variant: 'url', url: imageUrl });
      const result = await analyzeImage(imageUrl);
      if (onAnalysisComplete) {
        onAnalysisComplete(result);
      }
      setImageUrl('');
      } catch (error) {
        logger.error('URL structuralization failed:', error);
        setUrlError('Failed to structuralize from URL');
    } finally {
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

  const handleDrop = async (e) => {
    e.preventDefault();
    setIsDragOver(false);

    const dt = e.dataTransfer;
    if (!dt) return;

    // Prefer files first
    if (dt.files && dt.files.length > 0) {
      const file = dt.files[0];
      if (file && file.type && file.type.startsWith('image/')) {
        if (checkPaymentRequired()) return;
        setIsProcessing(true);
        setUrlError('');
        try {
          const reader = new FileReader();
          reader.onload = async (ev) => {
            try {
              const dataUrl = ev.target.result;
              logger.debug('Drop analyze payload summary', { variant: 'file-drop' });
              const result = await analyzeImage(dataUrl);
              if (onAnalysisComplete) onAnalysisComplete(result);
              setImageUrl('');
            } catch (err) {
              logger.error('File drop structuralization failed:', err);
              setUrlError('Failed to structuralize dropped image');
            } finally {
              setIsProcessing(false);
            }
          };
          reader.readAsDataURL(file);
        } catch (_) {
          setIsProcessing(false);
        }
      }
      return;
    }

    // Otherwise, attempt URL/text drop
    let dropped = dt.getData('text/uri-list') || dt.getData('text/plain');
    if (dropped && typeof dropped === 'string') {
      const url = dropped.trim();
      if (!url) return;
      setImageUrl(url);

      const validation = validateUrl(url);
      if (validation) {
        setUrlError(validation);
        return;
      }

      if (checkPaymentRequired()) return;

      setIsProcessing(true);
      setUrlError('');
      try {
        const result = await analyzeImage(url);
        if (onAnalysisComplete) onAnalysisComplete(result);
        setImageUrl('');
      } catch (error) {
        logger.error('URL drop structuralization failed:', error);
        setUrlError('Failed to structuralize from dropped URL');
      } finally {
        setIsProcessing(false);
      }
    }
  };

  return (
    <form onSubmit={handleSubmit}>
      <div 
        className={`url-row${isDragOver ? ' dragover' : ''}`}
        onDragOver={handleDragOver}
        onDragLeave={handleDragLeave}
        onDrop={handleDrop}
      >
        <input
          type="url"
          placeholder="Enter or drag image URL..."
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
          } else if (!molecularData?.skipCid && molecularData?.cid) {
            try {
              const cid = encodeURIComponent(String(molecularData.cid));
              const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/SDF?record_type=3d`;
              const fetched = await fetch(url);
              if (fetched.ok) {
                const text = await fetched.text();
                viewer.addModel(text, 'sdf');
              }
            } catch (_) {}
          } else if (molecularData.smiles) {
            try { viewer.addModel(molecularData.smiles, 'smiles'); }
            catch (_) { /* fall through to error path if load fails */ }
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
          // Fit, then back the camera off slightly to avoid being inside
          // very small/degenerate bounding boxes (e.g., single-atom SDFs).
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

const MolecularColumn = ({ column, onRemove, showRemove = true }) => {
  // Render column header and viewers
  
  return (
    <div className="column">
      {/* GUARANTEED VISIBLE HEADER */}
      <div className="column-header">
        <div className="column-meta">
          <div className="column-title">
            {column.query}
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

      {false && column.loading && column.viewers.length === 0 && (
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
  // Visual tests only run when explicitly enabled (not automatic in dev mode)
  // Enable via: REACT_APP_RUN_VISUAL_TESTS=true or URL ?visual-tests
  const [autoVisualMode, setAutoVisualMode] = useState(
    process.env.NODE_ENV === 'development' && 
    (process.env.REACT_APP_RUN_VISUAL_TESTS === 'true' || new URLSearchParams(window.location.search).has('visual-tests'))
  );
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

  // Handle keyboard shortcuts (desktop only)
  useEffect(() => {
    // Skip keyboard shortcuts on mobile devices
    if (isMobileDevice()) {
      logger.info('Skipping shortcuts - mobile device detected');
      return;
    }

    // Create keyboard handler with mode actions
    const actions = {
      focusInput: () => {
        const inputElement = document.getElementById('object-input');
        if (inputElement) {
          inputElement.focus();
          inputElement.select();
        }
      },
      cameraMode: () => {
        setCameraMode(true);
        setPhotoMode(false);
        setLinkMode(false);
      },
      photoMode: () => {
        setCameraMode(false);
        setPhotoMode(true);
        setLinkMode(false);
      },
      linkMode: () => {
        setCameraMode(false);
        setPhotoMode(false);
        setLinkMode(true);
      }
    };

    const keyboardHandler = createKeyboardHandler(actions);
    document.addEventListener('keydown', keyboardHandler, { capture: true });

    return () => document.removeEventListener('keydown', keyboardHandler, { capture: true });
  }, [setCameraMode, setPhotoMode, setLinkMode]);

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
              logger.error('SDF generation failed:', sdfError);
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
      logger.error('Analysis failed:', error);
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
          logger.error('SDF generation failed:', sdfError);
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
            logger.debug('Settings button clicked', { showSettings });
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
                  logger.info('Column mode changed', { value: e.target.value });
                  setColumnMode(e.target.value);
                }}
                className="select"
              >
                <option value="replace">
                  Replace (Default) - New analysis displaces previous
                </option>
                <option value="accumulate">
                  Accumulate - Add columns side by side
                </option>
              </select>
            </div>

            {process.env.NODE_ENV === 'development' && (
              <div className="settings-field">
                <label className="label">
                  <input
                    type="checkbox"
                    checked={autoVisualMode}
                    onChange={(e) => setAutoVisualMode(e.target.checked)}
                    className="checkbox-input"
                  />
                  Enable Visual Tests (Dev Mode)
                </label>
              </div>
            )}

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