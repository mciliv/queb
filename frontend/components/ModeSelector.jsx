import React from 'react';

const ModeSelector = ({ cameraMode, setCameraMode, photoMode, setPhotoMode }) => {
  const handleCameraChange = (e) => {
    setCameraMode(e.target.checked);
    if (e.target.checked) {
      setPhotoMode(false);
    }
  };

  const handlePhotoChange = (e) => {
    setPhotoMode(e.target.checked);
    if (e.target.checked) {
      setCameraMode(false);
    }
  };

  return (
    <div className="input-mode-section">
      <div className="mode-selector">
        <input 
          type="checkbox" 
          id="camera-mode" 
          checked={cameraMode}
          onChange={handleCameraChange}
        />
        <label htmlFor="camera-mode" className="mode-label">
          <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
            <circle cx="12" cy="12" r="5" fill="currentColor" opacity="0.8"/>
            <circle cx="12" cy="12" r="9" stroke="currentColor" fill="none"/>
          </svg>
          <span>Camera</span>
        </label>

        <input 
          type="checkbox" 
          id="photo-mode" 
          checked={photoMode}
          onChange={handlePhotoChange}
        />
        <label htmlFor="photo-mode" className="mode-label">
          <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
            <rect x="3" y="3" width="18" height="18" rx="2" ry="2"/>
            <circle cx="9" cy="9" r="2"/>
            <path d="M21 15l-3.086-3.086a2 2 0 0 0-2.828 0L6 21"/>
          </svg>
          <span>Photo</span>
        </label>
      </div>
    </div>
  );
};

export default ModeSelector;