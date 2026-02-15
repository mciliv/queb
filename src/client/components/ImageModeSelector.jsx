import React, { memo } from 'react';
import { isMobileDevice } from '../utils/device.js';

const cameraSvg = `
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
  <circle cx="12" cy="12" r="5" fill="currentColor" opacity="0.8" />
  <circle cx="12" cy="12" r="9" fill="none" stroke="currentColor" />
</svg>
`;

const photoSvg = `
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
  <rect x="3" y="3" width="18" height="18" rx="2" ry="2" />
  <circle cx="9" cy="9" r="2" />
  <path d="M21 15l-3.086-3.086a2 2 0 0 0-2.828 0L6 21" />
</svg>
`;

const linkSvg = `
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
  <path d="M10 13a5 5 0 0 0 7.54.54l3-3a5 5 0 0 0-7.07-7.07l-1.72 1.71" />
  <path d="M14 11a5 5 0 0 0-7.54-.54l-3 3a5 5 0 0 0 7.07 7.07l1.71-1.71" />
</svg>
`;

const ModeButton = memo(({ id, svgMarkup, label, shortcut, isActive, onClick }) => {
  const isMobile = isMobileDevice();

  return (
    <button
      className={`mode-btn${isActive ? ' active' : ''}`}
      onClick={onClick}
      title={isMobile ? label : `${label} (${shortcut})`}
    >
      <span
        className="mode-icon"
        aria-hidden="true"
        dangerouslySetInnerHTML={{ __html: svgMarkup }}
      />
      {!isMobile && <span className="mode-btn-shortcut">{shortcut}</span>}
    </button>
  );
});

ModeButton.displayName = 'ModeButton';

const ImageModeSelector = ({ cameraMode, setCameraMode, photoMode, setPhotoMode, linkMode, setLinkMode }) => {
  const handleModeSelect = (mode) => {
    setCameraMode(mode === 'camera');
    setPhotoMode(mode === 'photo');
    setLinkMode(mode === 'link');
  };

  return (
    <div className="mode-row">
      <ModeButton 
        id="camera" 
        svgMarkup={cameraSvg} 
        label="Capture from camera" 
        shortcut="Alt+1"
        isActive={cameraMode}
        onClick={() => handleModeSelect('camera')}
      />
      <ModeButton 
        id="photo" 
        svgMarkup={photoSvg} 
        label="Upload image" 
        shortcut="Alt+2"
        isActive={photoMode}
        onClick={() => handleModeSelect('photo')}
      />
      <ModeButton 
        id="link" 
        svgMarkup={linkSvg} 
        label="Enter image link" 
        shortcut="Alt+3"
        isActive={linkMode}
        onClick={() => handleModeSelect('link')}
      />
    </div>
  );
};

export default memo(ImageModeSelector);

