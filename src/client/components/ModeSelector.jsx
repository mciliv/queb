import React, { memo } from 'react';
import { useAppController } from '../AppController';
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

const ModeButton = memo(({ id, svgMarkup, label, shortcut }) => {
  const { activeMode, switchMode } = useAppController();
  const isMobile = isMobileDevice();

  return (
    <button
      className={`mode-btn${activeMode === id ? ' active' : ''}`}
      onClick={() => switchMode(id)}
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

const ModeSelector = () => (
  <div className="mode-row">
    <ModeButton id="camera" svgMarkup={cameraSvg} label="Capture from camera" shortcut="Alt+1" />
    <ModeButton id="photo" svgMarkup={photoSvg} label="Upload image" shortcut="Alt+2" />
    <ModeButton id="link" svgMarkup={linkSvg} label="Enter image link" shortcut="Alt+3" />
  </div>
);

export default memo(ModeSelector);


