import React, { useState, useEffect, useRef, useCallback, useMemo } from 'react';
import { VALIDATION_PATTERNS } from '../utils/config-loader.js';
import { isMobileDevice } from '../utils/device.js';

const TextInput = ({ value, onChange, onSubmit, isProcessing, error }) => {
  const [localError, setLocalError] = useState('');
  const [isValidating, setIsValidating] = useState(false);
  const [suggestions, setSuggestions] = useState([]);
  const [showSuggestions, setShowSuggestions] = useState(false);
  const [selectedSuggestionIndex, setSelectedSuggestionIndex] = useState(-1);
  const inputRef = useRef(null);
  const lastTriggerTimeRef = useRef(0);
  const debounceTimeoutRef = useRef(null);
  const isMac = navigator.userAgent.toLowerCase().includes('mac');
  const keyboardHint = isMac ? '⌘K' : 'Ctrl+K';

  // Common molecule suggestions for autocomplete
  const moleculeSuggestions = useMemo(() => [
    'caffeine', 'glucose', 'water', 'ethanol', 'acetaminophen', 'aspirin',
    'vitamin C', 'vitamin D', 'insulin', 'dopamine', 'serotonin', 'adrenaline',
    'cholesterol', 'testosterone', 'estrogen', 'THC', 'nicotine', 'morphine',
    'penicillin', 'ibuprofen', 'sodium chloride', 'carbon dioxide', 'methane',
    'benzene', 'toluene', 'acetone', 'formaldehyde', 'ammonia', 'hydrogen peroxide'
  ], []);

  useEffect(() => {
    if (localError && value) {
      setLocalError('');
    }
  }, [value, localError]);

  // Debounced input change handler
  const debouncedOnChange = useCallback((newValue) => {
    // Clear previous timeout
    if (debounceTimeoutRef.current) {
      clearTimeout(debounceTimeoutRef.current);
    }

    // Immediate update for UI responsiveness
    onChange(newValue);

    // Debounced suggestions update
    debounceTimeoutRef.current = setTimeout(() => {
      if (newValue.trim().length >= 2) {
        const filteredSuggestions = moleculeSuggestions
          .filter(suggestion => 
            suggestion.toLowerCase().includes(newValue.toLowerCase()) &&
            suggestion.toLowerCase() !== newValue.toLowerCase()
          )
          .slice(0, 5);
        setSuggestions(filteredSuggestions);
        setShowSuggestions(filteredSuggestions.length > 0);
        setSelectedSuggestionIndex(-1);
      } else {
        setShowSuggestions(false);
        setSuggestions([]);
      }
    }, 150); // 150ms debounce
  }, [onChange, moleculeSuggestions]);

  // Clean up debounce timeout on unmount
  useEffect(() => {
    return () => {
      if (debounceTimeoutRef.current) {
        clearTimeout(debounceTimeoutRef.current);
      }
    };
  }, []);

  // Handle suggestion selection
  const selectSuggestion = useCallback((suggestion) => {
    onChange(suggestion);
    setShowSuggestions(false);
    setSuggestions([]);
    setSelectedSuggestionIndex(-1);
    if (inputRef.current) {
      inputRef.current.focus();
    }
  }, [onChange]);

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
    if (showSuggestions && suggestions.length > 0) {
      if (e.key === 'ArrowDown') {
        e.preventDefault();
        setSelectedSuggestionIndex(prev => 
          prev < suggestions.length - 1 ? prev + 1 : 0
        );
        return;
      }
      if (e.key === 'ArrowUp') {
        e.preventDefault();
        setSelectedSuggestionIndex(prev => 
          prev > 0 ? prev - 1 : suggestions.length - 1
        );
        return;
      }
      if (e.key === 'Tab' && selectedSuggestionIndex >= 0) {
        e.preventDefault();
        selectSuggestion(suggestions[selectedSuggestionIndex]);
        return;
      }
      if (e.key === 'Escape') {
        setShowSuggestions(false);
        setSelectedSuggestionIndex(-1);
        return;
      }
    }
    
    if (e.key === 'Enter') {
      e.preventDefault();
      e.stopPropagation();
      if (showSuggestions && selectedSuggestionIndex >= 0) {
        selectSuggestion(suggestions[selectedSuggestionIndex]);
      } else {
        await handleSubmit();
      }
    }
  };

  // Handle input blur to hide suggestions
  const handleBlur = useCallback((e) => {
    // Delay hiding to allow clicking on suggestions
    setTimeout(() => {
      if (!e.currentTarget.contains(document.activeElement)) {
        setShowSuggestions(false);
        setSelectedSuggestionIndex(-1);
      }
    }, 150);
  }, []);

  const displayError = localError || error;

  return (
    <div className="input-wrapper">
      <div className="input-row">
        <div className="input-container" style={{ position: 'relative' }}>
          <input
            ref={inputRef}
            id="object-input"
            type="text"
            placeholder="Specify object... (start typing for suggestions)"
            className={`input-base${displayError ? ' input-error' : ''}`}
            value={value}
            onChange={(e) => debouncedOnChange(e.target.value)}
            onKeyDown={handleKeyDown}
            onBlur={handleBlur}
            aria-describedby={displayError ? 'input-error' : undefined}
            aria-expanded={showSuggestions}
            aria-haspopup="listbox"
            autoComplete="off"
          />
          
          {showSuggestions && suggestions.length > 0 && (
            <div 
              className="suggestions-dropdown"
              style={{
                position: 'absolute',
                top: '100%',
                left: 0,
                right: 0,
                backgroundColor: 'white',
                border: '1px solid #ddd',
                borderTop: 'none',
                borderRadius: '0 0 4px 4px',
                boxShadow: '0 2px 8px rgba(0,0,0,0.1)',
                zIndex: 1000,
                maxHeight: '200px',
                overflowY: 'auto'
              }}
              role="listbox"
            >
              {suggestions.map((suggestion, index) => (
                <div
                  key={suggestion}
                  className={`suggestion-item ${index === selectedSuggestionIndex ? 'selected' : ''}`}
                  style={{
                    padding: '8px 12px',
                    cursor: 'pointer',
                    backgroundColor: index === selectedSuggestionIndex ? '#f0f0f0' : 'transparent',
                    borderBottom: index < suggestions.length - 1 ? '1px solid #eee' : 'none'
                  }}
                  onClick={() => selectSuggestion(suggestion)}
                  onMouseEnter={() => setSelectedSuggestionIndex(index)}
                  role="option"
                  aria-selected={index === selectedSuggestionIndex}
                >
                  {suggestion}
                </div>
              ))}
            </div>
          )}
        </div>
        
        {value.trim() && (
          <button 
            className="btn-icon"
            onClick={handleSubmit}
            aria-label="Structuralize"
          >
            →
         </button>
        )}
      </div>
      
      {!value.trim() && !isMobileDevice() && (
        <div className="kbd-hint">{keyboardHint}</div>
      )}
      
      {displayError && (
        <div id="input-error" className="error-text" role="alert">
          {displayError}
        </div>
      )}
    </div>
  );
};

export default TextInput;