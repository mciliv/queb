import React, { useRef, useEffect, useState } from 'react';
import { usePayment } from './PaymentContext';
import { useApi } from '../hooks/useApi';

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
    <div className="camera-container" onClick={handleCameraClick}>
      <video 
        ref={videoRef}
        autoPlay 
        playsInline 
        muted 
        style={{ width: '100%' }}
      />
      
      {!hasPermission && (
        <div className="permission-message">{permissionMessage}</div>
      )}
      
      <div className="crosshair"></div>
      
      {showSwitchCamera && (
        <div className="switch-camera-container">
          <button 
            type="button" 
            className="switch-camera-btn"
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
      
      <div className="instruction-text">
        Center object in circle & tap, or type name above
      </div>
    </div>
  );
};

export default CameraSection;