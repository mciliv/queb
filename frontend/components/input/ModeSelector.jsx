import React from 'react';

const ModeSelector = ({ cameraMode, setCameraMode, photoMode, setPhotoMode, linkMode, setLinkMode }) => {
  const handleModeSelect = (mode) => {
    // Reset all modes
    setCameraMode(false);
    setPhotoMode(false);
    if (setLinkMode) setLinkMode(false);
    
    // Set selected mode
    switch(mode) {
      case 'camera':
        setCameraMode(true);
        break;
      case 'photo':
        setPhotoMode(true);
        break;
      case 'link':
        if (setLinkMode) setLinkMode(true);
        break;
    }
  };

  const styles = {
    container: {
      display: 'flex',
      gap: '10px',
      marginTop: '10px'
    },
    button: {
      background: 'transparent',
      border: '1px solid rgba(255, 255, 255, 0.2)',
      color: '#ffffff',
      padding: '8px 16px',
      borderRadius: '4px',
      cursor: 'pointer',
      display: 'flex',
      alignItems: 'center',
      gap: '8px',
      fontSize: '13px',
      transition: 'all 0.2s'
    },
    activeButton: {
      background: 'rgba(255, 255, 255, 0.1)',
      borderColor: 'rgba(255, 255, 255, 0.4)'
    }
  };

  return (
    <div style={styles.container}>
      <button 
        style={{...styles.button, ...(cameraMode ? styles.activeButton : {})}}
        onClick={() => handleModeSelect('camera')}
        title="Capture from camera"
      >
        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
          <circle cx="12" cy="12" r="5" fill="currentColor" opacity="0.8"/>
          <circle cx="12" cy="12" r="9" stroke="currentColor" fill="none"/>
        </svg>
        Camera
      </button>

      <button 
        style={{...styles.button, ...(photoMode ? styles.activeButton : {})}}
        onClick={() => handleModeSelect('photo')}
        title="Upload image"
      >
        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
          <rect x="3" y="3" width="18" height="18" rx="2" ry="2"/>
          <circle cx="9" cy="9" r="2"/>
          <path d="M21 15l-3.086-3.086a2 2 0 0 0-2.828 0L6 21"/>
        </svg>
        Image
      </button>

      <button 
        style={{...styles.button, ...(linkMode ? styles.activeButton : {})}}
        onClick={() => handleModeSelect('link')}
        title="Enter image link"
      >
        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
          <path d="M10 13a5 5 0 0 0 7.54.54l3-3a5 5 0 0 0-7.07-7.07l-1.72 1.71"/>
          <path d="M14 11a5 5 0 0 0-7.54-.54l-3 3a5 5 0 0 0 7.07 7.07l1.71-1.71"/>
        </svg>
        Link
      </button>
    </div>
  );
};

export default ModeSelector;