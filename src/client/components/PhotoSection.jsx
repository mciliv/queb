import React, { useState, useRef, useEffect } from 'react';
import { usePayment } from './ui/PaymentContext';
import { useApi } from '../hooks/useApi';
import { isMobileDevice } from '../utils/device.js';

const PhotoSection = ({ isProcessing, setIsProcessing, setCurrentAnalysisType, onAnalysisComplete, openOnMount = false }) => {
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

              const result = await analyzeImage(dataUrl);
              if (onAnalysisComplete) {
                onAnalysisComplete(result);
              }
            } catch (innerErr) {
              console.error('Photo header crop failed:', innerErr);
            } finally {
              setIsProcessing(false);
            }
          };
          img.onerror = () => {
            // Fallback: pass original if crop fails to load
            analyzeImage(dataUrl)
              .then((result) => {
                if (onAnalysisComplete) {
                  onAnalysisComplete(result);
                }
              })
              .catch((err) => console.error('Photo structuralization failed:', err))
              .finally(() => setIsProcessing(false));
          };
          img.src = dataUrl;
        } catch (error) {
          console.error('File reading or processing failed:', error);
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

  // If openOnMount is set, immediately open file chooser when section appears
  useEffect(() => {
    if (openOnMount && fileInputRef.current) {
      fileInputRef.current.click();
    }
  }, [openOnMount]);


  // Mobile: Simple button interface
  if (isMobile) {
    return (
      <div>
        <input
          ref={fileInputRef}
          type="file"
          accept="image/*"
          capture="environment"
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
        capture="environment"
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

export default PhotoSection;
