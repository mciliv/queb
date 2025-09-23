import React from 'react';
import { isMobileDevice } from '../utils/device.js';

const ModeSelector = ({ mode, setMode }) => {
  const isMobile = isMobileDevice();
  const handleModeSelect = (next) => setMode(next);

  return (
    <div className="mode-row">
      <button 
        className={`mode-btn${mode === 'camera' ? ' active' : ''}`}
        onClick={() => handleModeSelect('camera')}
        title={isMobile ? "Capture from camera" : "Capture from camera (⌘⇧1)"}
      >
        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
          <circle cx="12" cy="12" r="5" fill="currentColor" opacity="0.8"/>
          <circle cx="12" cy="12" r="9" stroke="currentColor" fill="none"/>
        </svg>
        {!isMobile && (
          <span className="mode-btn-shortcut">
            {navigator.userAgent.toUpperCase().indexOf('MAC') >= 0 ? '⌘⇧1' : 'Ctrl+Shift+1'}
          </span>
        )}
      </button>

      <button 
        className={`mode-btn${mode === 'photo' ? ' active' : ''}`}
        onClick={() => handleModeSelect('photo')}
        title={isMobile ? "Upload image" : "Upload image (⌘⇧2)"}
      >
        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
          <rect x="3" y="3" width="18" height="18" rx="2" ry="2"/>
          <circle cx="9" cy="9" r="2"/>
          <path d="M21 15l-3.086-3.086a2 2 0 0 0-2.828 0L6 21"/>
        </svg>
        {!isMobile && (
          <span className="mode-btn-shortcut">
            {navigator.userAgent.toUpperCase().indexOf('MAC') >= 0 ? '⌘⇧2' : 'Ctrl+Shift+2'}
          </span>
        )}
      </button>

      <button 
        className={`mode-btn${mode === 'link' ? ' active' : ''}`}
        onClick={() => handleModeSelect('link')}
        title={isMobile ? "Enter image link" : "Enter image link (⌘⇧3)"}
      >
        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
          <path d="M10 13a5 5 0 0 0 7.54.54l3-3a5 5 0 0 0-7.07-7.07l-1.72 1.71"/>
          <path d="M14 11a5 5 0 0 0-7.54-.54l-3 3a5 5 0 0 0 7.07 7.07l1.71-1.71"/>
        </svg>
        {!isMobile && (
          <span className="mode-btn-shortcut">
            {navigator.userAgent.toUpperCase().indexOf('MAC') >= 0 ? '⌘⇧3' : 'Ctrl+Shift+3'}
          </span>
        )}
      </button>
    </div>
  );
};

export default ModeSelector;
