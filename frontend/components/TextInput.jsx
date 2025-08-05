import React, { useState, useEffect, useRef } from 'react';

const TextInput = ({ value, onChange, onSubmit, isProcessing, error }) => {
  const [localError, setLocalError] = useState('');
  const [isValidating, setIsValidating] = useState(false);
  const inputRef = useRef(null);
  const isMac = navigator.platform.toLowerCase().includes('mac');
  const keyboardHint = isMac ? '⌘K' : 'Ctrl+K';

  // Clear local error when value changes
  useEffect(() => {
    if (localError && value) {
      setLocalError('');
    }
  }, [value, localError]);

  // Focus input when component mounts
  useEffect(() => {
    if (inputRef.current) {
      inputRef.current.focus();
    }
  }, []);

  const validateInput = (text) => {
    if (!text || !text.trim()) {
      return 'Please enter a molecule or object to analyze';
    }
    if (text.trim().length < 2) {
      return 'Input must be at least 2 characters long';
    }
    if (text.trim().length > 500) {
      return 'Input must be less than 500 characters';
    }
    return null;
  };

  const handleKeyDown = async (e) => {
    if (e.key === 'Enter' && !isProcessing) {
      e.preventDefault();
      e.stopPropagation();
      
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
        setLocalError(err.message || 'Analysis failed. Please try again.');
      } finally {
        setIsValidating(false);
      }
    }
  };

  const handleSubmit = async () => {
    if (isProcessing || isValidating) return;
    
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
      setLocalError(err.message || 'Analysis failed. Please try again.');
    } finally {
      setIsValidating(false);
    }
  };

  const displayError = localError || error;
  const isDisabled = isProcessing || isValidating;

  return (
    <div className="top-bar">
      <div className="input-container">
        <input
          ref={inputRef}
          id="object-input"
          type="text"
          placeholder="Specify an object or molecule (e.g., 'water', 'aspirin', 'CCO')..."
          className={`text-input ${displayError ? 'error' : ''}`}
          value={value}
          onChange={(e) => onChange(e.target.value)}
          onKeyDown={handleKeyDown}
          disabled={isDisabled}
          aria-describedby={displayError ? 'input-error' : undefined}
        />
        <div className="keyboard-hint">
          <span className="hint-key">{keyboardHint}</span>
        </div>
        {!isDisabled && value.trim() && (
          <button 
            className="submit-button"
            onClick={handleSubmit}
            aria-label="Submit analysis"
          >
            →
          </button>
        )}
      </div>
      {displayError && (
        <div id="input-error" className="input-error" role="alert">
          {displayError}
        </div>
      )}
      {(isProcessing || isValidating) && (
        <div className="processing-indicator">
          <span className="spinner"></span>
          Analyzing...
        </div>
      )}
    </div>
  );
};

export default TextInput;