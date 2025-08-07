import React, { useRef, useEffect, useState } from 'react';
import { usePayment } from '../ui/PaymentContext';
import { useApi } from '../../hooks/useApi';

// Inline styles for camera component
const styles = {
  cameraContainer: {
    position: 'relative',
    width: '100%',
    maxWidth: '400px',
    borderRadius: '8px',
    overflow: 'hidden',
    cursor: 'pointer',
    background: '#000',
    aspectRatio: '4/3'
  },
  permissionMessage: {
    position: 'absolute',
    top: '50%',
    left: '50%',
    transform: 'translate(-50%, -50%)',
    color: '#ffffff',
    textAlign: 'center',
    padding: '20px',
    background: 'rgba(0, 0, 0, 0.8)',
    borderRadius: '8px',
    fontSize: '14px'
  },
  crosshair: {
    position: 'absolute',
    top: '50%',
    left: '50%',
    transform: 'translate(-50%, -50%)',
    width: '60px',
    height: '60px',
    border: '2px solid rgba(255, 255, 255, 0.8)',
    borderRadius: '50%',
    pointerEvents: 'none'
  },
  switchCameraContainer: {
    position: 'absolute',
    top: '16px',
    right: '16px'
  },
  switchCameraBtn: {
    display: 'flex',
    alignItems: 'center',
    gap: '8px',
    background: 'rgba(0, 0, 0, 0.7)',
    color: '#ffffff',
    border: 'none',
    borderRadius: '20px',
    padding: '8px 12px',
    fontSize: '12px',
    cursor: 'pointer'
  },
  instructionText: {
    position: 'absolute',
    bottom: '16px',
    left: '50%',
    transform: 'translateX(-50%)',
    color: '#ffffff',
    fontSize: '12px',
    textAlign: 'center',
    background: 'rgba(0, 0, 0, 0.7)',
    padding: '8px 12px',
    borderRadius: '16px',
    whiteSpace: 'nowrap'
  }
};

const CameraSection = ({ isProcessing, setIsProcessing, setCurrentAnalysisType, onAnalysisComplete }) => {
  const videoRef = useRef(null);
  const canvasRef = useRef(null);
  const [hasPermission, setHasPermission] = useState(false);
  const [permissionMessage, setPermissionMessage] = useState('');
  const [showSwitchCamera, setShowSwitchCamera] = useState(false);
  const [stream, setStream] = useState(null);
  const { checkPaymentRequired } = usePayment();
  const { analyzeImage } = useApi();

  useEffect(() => {
    requestCameraAccess();
    return () => {
      if (stream) {
        stream.getTracks().forEach(track => track.stop());
      }
    };
  }, []);

  const requestCameraAccess = async () => {
    try {
      const mediaStream = await navigator.mediaDevices.getUserMedia({ 
        video: { facingMode: 'environment' },
        audio: false 
      });
      
      if (videoRef.current) {
        videoRef.current.srcObject = mediaStream;
      }
      
      setStream(mediaStream);
      setHasPermission(true);
      
      // Check if multiple cameras available
      const devices = await navigator.mediaDevices.enumerateDevices();
      const videoDevices = devices.filter(device => device.kind === 'videoinput');
      setShowSwitchCamera(videoDevices.length > 1);
    } catch (error) {
      console.error('Camera access error:', error);
      setPermissionMessage('Camera access required for molecular analysis');
      setHasPermission(false);
    }
  };

  const handleCameraClick = async () => {
    if (isProcessing || !hasPermission || !videoRef.current) return;

    if (checkPaymentRequired()) {
      return;
    }

    setIsProcessing(true);
    setCurrentAnalysisType('camera');
    
    try {
      // Create canvas and capture frame
      const canvas = document.createElement('canvas');
      const video = videoRef.current;
      canvas.width = video.videoWidth;
      canvas.height = video.videoHeight;
      const ctx = canvas.getContext('2d');
      ctx.drawImage(video, 0, 0);
      
      // Convert to base64
      const imageData = canvas.toDataURL('image/jpeg', 0.8);
      
      // Analyze the image
      const result = await analyzeImage(imageData, 'Camera capture');
      
      if (onAnalysisComplete) {
        onAnalysisComplete(result);
      }
    } catch (error) {
      console.error('Camera analysis failed:', error);
    } finally {
      setIsProcessing(false);
    }
  };

  const switchCamera = async () => {
    if (stream) {
      stream.getTracks().forEach(track => track.stop());
    }
    
    // Toggle between front and back camera
    const currentFacingMode = stream?.getVideoTracks()[0]?.getSettings()?.facingMode || 'environment';
    const newFacingMode = currentFacingMode === 'environment' ? 'user' : 'environment';
    
    try {
      const mediaStream = await navigator.mediaDevices.getUserMedia({ 
        video: { facingMode: newFacingMode },
        audio: false 
      });
      
      if (videoRef.current) {
        videoRef.current.srcObject = mediaStream;
      }
      
      setStream(mediaStream);
    } catch (error) {
      console.error('Camera switch error:', error);
    }
  };

  return (
    <div style={styles.cameraContainer} onClick={handleCameraClick}>
      <video 
        ref={videoRef}
        autoPlay 
        playsInline 
        muted 
        style={{ width: '100%', height: '100%', objectFit: 'cover' }}
      />
      
      {!hasPermission && (
        <div style={styles.permissionMessage}>{permissionMessage}</div>
      )}
      
      <div style={styles.crosshair}></div>
      
      {showSwitchCamera && (
        <div style={styles.switchCameraContainer}>
          <button 
            type="button" 
            style={styles.switchCameraBtn}
            onClick={(e) => {
              e.stopPropagation();
              switchCamera();
            }}
          >
            <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
              <circle cx="12" cy="12" r="5" fill="currentColor" opacity="0.8"/>
              <circle cx="12" cy="12" r="9" stroke="currentColor" fill="none"/>
            </svg>
            <span>Switch Camera</span>
          </button>
        </div>
      )}
      
      <div style={styles.instructionText}>
        Center object in circle & tap, or type name above
      </div>
    </div>
  );
};

export default CameraSection;