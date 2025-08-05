import React, { useEffect } from 'react';

const TextInput = ({ value, onChange, onSubmit, isProcessing }) => {
  const isMac = navigator.platform.toLowerCase().includes('mac');
  const keyboardHint = isMac ? '⌘K' : 'Ctrl+K';

  const handleKeyDown = (e) => {
    if (e.key === 'Enter' && !isProcessing && value.trim()) {
      e.preventDefault();
      e.stopPropagation();
      onSubmit(value);
    }
  };

  return (
    <div className="top-bar">
      <div className="input-container">
        <input
          id="object-input"
          type="text"
          placeholder="Specify an object..."
          className="text-input"
          value={value}
          onChange={(e) => onChange(e.target.value)}
          onKeyDown={handleKeyDown}
          disabled={isProcessing}
        />
        <div className="keyboard-hint">
          <span className="hint-key">{keyboardHint}</span>
        </div>
      </div>
    </div>
  );
};

export default TextInput;