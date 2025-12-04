import React, { useState, useEffect, useRef } from 'react';
import { VALIDATION_PATTERNS } from '../utils/config-loader.js';
import { isMobileDevice } from '../utils/device.js';

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
          placeholder="What's in..."
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

export default TextInput;
