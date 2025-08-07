import React, { useState, useEffect, useRef } from 'react';

// Inline styles for text input component
const styles = {
  topBar: {
    position: 'relative',
    width: '100%',
    marginBottom: '20px'
  },
  inputContainer: {
    position: 'relative',
    display: 'flex',
    alignItems: 'center',
    gap: '8px'
  },
  textInput: {
    width: '100%',
    padding: '12px 16px',
    background: 'rgba(255, 255, 255, 0.08)',
    border: 'none',
    borderRadius: '8px',
    color: '#ffffff',
    fontSize: '14px',
    outline: 'none',
    fontFamily: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif'
  },
  textInputError: {
    width: '100%',
    padding: '12px 16px',
    background: 'rgba(255, 255, 255, 0.08)',
    border: '1px solid rgba(255, 100, 100, 0.5)',
    borderRadius: '8px',
    color: '#ffffff',
    fontSize: '14px',
    outline: 'none',
    fontFamily: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif'
  },
  keyboardHint: {
    position: 'absolute',
    right: '50px',
    top: '50%',
    transform: 'translateY(-50%)',
    pointerEvents: 'none'
  },
  hintKey: {
    background: 'rgba(255, 255, 255, 0.1)',
    color: 'rgba(255, 255, 255, 0.6)',
    padding: '2px 6px',
    borderRadius: '4px',
    fontSize: '11px',
    fontWeight: '500'
  },
  submitButton: {
    background: 'rgba(255, 255, 255, 0.1)',
    border: 'none',
    borderRadius: '6px',
    color: '#ffffff',
    padding: '8px 12px',
    cursor: 'pointer',
    fontSize: '16px',
    minWidth: '40px',
    height: '40px',
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center'
  },
  inputError: {
    color: '#ff6b6b',
    fontSize: '12px',
    marginTop: '8px',
    padding: '0 4px'
  },
  processingIndicator: {
    display: 'flex',
    alignItems: 'center',
    gap: '8px',
    color: 'rgba(255, 255, 255, 0.8)',
    fontSize: '12px',
    marginTop: '8px',
    padding: '0 4px'
  },
  spinner: {
    width: '12px',
    height: '12px',
    border: '2px solid rgba(255, 255, 255, 0.3)',
    borderTop: '2px solid #ffffff',
    borderRadius: '50%',
    animation: 'spin 1s linear infinite'
  }
};

const TextInput = ({ value, onChange, onSubmit, isProcessing, error }) => {
  const [localError, setLocalError] = useState('');
  const [isValidating, setIsValidating] = useState(false);
  const inputRef = useRef(null);
  const lastTriggerTimeRef = useRef(0);
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
      
      // Prevent multiple rapid triggers with debouncing
      const now = Date.now();
      const DEBOUNCE_MS = 1000; // 1 second cooldown
      
      if (now - lastTriggerTimeRef.current < DEBOUNCE_MS) {
        console.log('🚫 Enter key debounced - too soon since last trigger');
        return;
      }
      
      lastTriggerTimeRef.current = now;
      
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
    <div style={styles.topBar}>
      <div style={styles.inputContainer}>
        <input
          ref={inputRef}
          id="object-input"
          type="text"
          placeholder="Specify an object or molecule (e.g., 'water', 'aspirin', 'CCO')..."
          style={displayError ? styles.textInputError : styles.textInput}
          value={value}
          onChange={(e) => onChange(e.target.value)}
          onKeyDown={handleKeyDown}
          disabled={isDisabled}
          aria-describedby={displayError ? 'input-error' : undefined}
        />
        <div style={styles.keyboardHint}>
          <span style={styles.hintKey}>{keyboardHint}</span>
        </div>
        {!isDisabled && value.trim() && (
          <button 
            style={styles.submitButton}
            onClick={handleSubmit}
            aria-label="Submit analysis"
          >
            →
          </button>
        )}
      </div>
      {displayError && (
        <div id="input-error" style={styles.inputError} role="alert">
          {displayError}
        </div>
      )}
      {(isProcessing || isValidating) && (
        <div style={styles.processingIndicator}>
          <span style={styles.spinner}></span>
          Analyzing...
        </div>
      )}
    </div>
  );
};

export default TextInput;