import React, { useState, useRef, useEffect } from 'react';
import { useApi } from '../hooks/useApi';

const LinkSection = ({ isProcessing, setIsProcessing, onPredictionComplete }) => {
  const [imageUrl, setImageUrl] = useState('');
  const [urlError, setUrlError] = useState('');
  const [isDragOver, setIsDragOver] = useState(false);
  const { analyzeImage } = useApi();
  const urlInputRef = useRef(null);

  // Auto-focus the URL input when component mounts
  useEffect(() => {
    if (urlInputRef.current) {
      urlInputRef.current.focus();
    }
  }, []);

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


    setIsProcessing(true);
    setUrlError('');
    
    try {
      const result = await analyzeImage(imageUrl);
      if (onPredictionComplete) {
        onPredictionComplete(result);
      }
      setImageUrl('');
      } catch (error) {
        console.error('URL structuralization failed:', error);
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
        setIsProcessing(true);
        setUrlError('');
        try {
          const reader = new FileReader();
          reader.onload = async (ev) => {
            try {
              const dataUrl = ev.target.result;
              const result = await analyzeImage(dataUrl);
              if (onPredictionComplete) onPredictionComplete(result);
              setImageUrl('');
            } catch (err) {
              console.error('File drop structuralization failed:', err);
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


      setIsProcessing(true);
      setUrlError('');
      try {
        const result = await analyzeImage(url);
        if (onPredictionComplete) onPredictionComplete(result);
        setImageUrl('');
      } catch (error) {
        console.error('URL drop structuralization failed:', error);
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
          ref={urlInputRef}
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

export default LinkSection;