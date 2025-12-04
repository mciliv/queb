/**
 * CameraSection.jsx - Live camera analysis component
 * 
 * This component enables real-time molecular analysis using the device camera.
 * Users can click on objects in the camera feed to analyze their chemical composition.
 * The component handles camera permissions, capture, and crop visualization.
 * 
 * Key features:
 * - Multi-camera support (standard + wide angle)
 * - Click-to-analyze with visual feedback
 * - Smart cropping around click point (25% of image)
 * - Fallback to uploaded photos if camera unavailable
 * - Responsive design for mobile and desktop
 */

import React, { useState, useEffect, useRef } from 'react';
import { useApi } from '../hooks/useApi';

/**
 * CameraSection Component
 * 
 * @param {Object} props
 * @param {boolean} props.isProcessing - Global processing state
 * @param {Function} props.setIsProcessing - Update processing state
 * @param {Function} props.setCurrentAnalysisType - Track analysis type for UI
 * @param {Function} props.onAnalysisComplete - Callback when analysis finishes
 */
const CameraSection = ({ isProcessing, setIsProcessing, setCurrentAnalysisType, onAnalysisComplete }) => {
  // Video element refs
  const videoRef = useRef(null);              // Main camera video element
  const wideVideoRef = useRef(null);          // Wide-angle camera (if available)
  
  // Camera state
  const [hasPermission, setHasPermission] = useState(false);
  const [permissionMessage, setPermissionMessage] = useState('');
  const [showSwitchCamera, setShowSwitchCamera] = useState(false);
  const [stream, setStream] = useState(null);           // MediaStream for main camera
  const [wideStream, setWideStream] = useState(null);   // MediaStream for wide camera
  const [availableCameras, setAvailableCameras] = useState([]);
  const [hasMultipleCameras, setHasMultipleCameras] = useState(false);
  
  // Click analysis state
  const [clickPosition, setClickPosition] = useState(null);  // Where user clicked
  const [showOutline, setShowOutline] = useState(false);     // Show crop outline
  const outlineTimeoutRef = useRef(null);                    // Timer for outline
  
  // External hooks
  const { analyzeImage } = useApi();

  useEffect(() => {
    requestCameraAccess();
    return () => {
      if (stream) {
        stream.getTracks().forEach(track => track.stop());
      }
      if (wideStream) {
        wideStream.getTracks().forEach(track => track.stop());
      }
      if (outlineTimeoutRef.current) {
        try { clearTimeout(outlineTimeoutRef.current); } catch (_) {}
      }
    };
  }, []);

  /**
   * Request camera permissions and initialize video streams
   * Attempts to get both standard and wide-angle cameras if available
   */
  const requestCameraAccess = async () => {
    try {
      // Enumerate all video input devices
      const devices = await navigator.mediaDevices.enumerateDevices();
      const videoDevices = devices.filter(device => device.kind === 'videoinput');
      setAvailableCameras(videoDevices);
      setShowSwitchCamera(videoDevices.length > 1);
      
      // Request standard camera with optimal settings
      const mediaStream = await navigator.mediaDevices.getUserMedia({ 
        video: { 
          facingMode: 'environment',      // Prefer rear camera on mobile
          width: { ideal: 1920 },         // Full HD for better analysis
          height: { ideal: 1080 },
          aspectRatio: { ideal: 16/9 }    // Standard widescreen
        },
        audio: false                      // No audio needed
      });
      
      if (videoRef.current) {
        videoRef.current.srcObject = mediaStream;
      }
      
      setStream(mediaStream);
      setHasPermission(true);
      
      // Try to get wide-angle camera if available
      try {
        const wideMediaStream = await navigator.mediaDevices.getUserMedia({
          video: {
            facingMode: 'environment',
            width: { ideal: 1920 },
            height: { ideal: 1080 },
            aspectRatio: { ideal: 16/9 },
            advanced: [{ zoom: { min: 0.5, max: 0.7 } }]
          },
          audio: false
        });
        
        if (wideVideoRef.current) {
          wideVideoRef.current.srcObject = wideMediaStream;
        }
        
        setWideStream(wideMediaStream);
        setHasMultipleCameras(true);
      } catch (wideError) {
        console.log('Wide camera not available:', wideError);
      }
      
    } catch (error) {
      console.error('Camera access error:', error);
      setPermissionMessage('Camera access required to structuralize from camera');
      setHasPermission(false);
    }
  };

  /**
   * Handle click on camera feed to analyze object at that location
   * Captures frame, calculates crop region, and sends for analysis
   * 
   * @param {MouseEvent|TouchEvent} e - Click/touch event with coordinates
   */
  const handleCameraClick = async (e) => {
    // Check permissions first
    if (!hasPermission) {
      await requestCameraAccess();
      return;
    }
    
    // Validate state
    if (!videoRef.current || isProcessing) return;

    setIsProcessing(true);
    setCurrentAnalysisType('camera');
    
    try {
      const canvas = document.createElement('canvas');
      const video = videoRef.current;
      const width = video.videoWidth;
      const height = video.videoHeight;
      canvas.width = width;
      canvas.height = height;
      const ctx = canvas.getContext('2d');
      
      // Capture wide camera if available
      let wideCanvas = null;
      if (hasMultipleCameras && wideVideoRef.current) {
        wideCanvas = document.createElement('canvas');
        const wideVideo = wideVideoRef.current;
        wideCanvas.width = wideVideo.videoWidth;
        wideCanvas.height = wideVideo.videoHeight;
        const wideCtx = wideCanvas.getContext('2d');
        wideCtx.drawImage(wideVideo, 0, 0);
      }

      // Determine click position relative to video frame (accounting for object-fit: cover)
      const container = e.currentTarget;
      const rect = container.getBoundingClientRect();
      const clickX = e.clientX - rect.left;
      const clickY = e.clientY - rect.top;

      const videoAspect = width / height;
      const containerAspect = rect.width / rect.height;
      let scale = 1;
      let offsetX = 0;
      let offsetY = 0;
      if (containerAspect > videoAspect) {
        // Container is wider; height overflows
        scale = rect.width / width;
        const scaledHeight = rect.width / videoAspect;
        offsetY = (rect.height - scaledHeight) / 2;
      } else {
        // Container is taller; width overflows
        scale = rect.height / height;
        const scaledWidth = rect.height * videoAspect;
        offsetX = (rect.width - scaledWidth) / 2;
      }

      // Map click to video coordinates
      const videoClickX = (clickX - offsetX) / scale;
      const videoClickY = (clickY - offsetY) / scale;

      // Crop square side (25% of min dimension)
      const cropSide = Math.floor(Math.min(width, height) * 0.25);
      const half = Math.floor(cropSide / 2);
      // Clamp center so crop stays within frame
      const centerX = Math.min(Math.max(videoClickX, half), width - half);
      const centerY = Math.min(Math.max(videoClickY, half), height - half);
      const cropX = Math.max(0, Math.floor(centerX - half));
      const cropY = Math.max(0, Math.floor(centerY - half));

      ctx.drawImage(video, 0, 0);
      const fullImageDataUrl = canvas.toDataURL('image/jpeg', 0.8);

      // Get wide image data if available
      const wideImageDataUrl = wideCanvas ? wideCanvas.toDataURL('image/jpeg', 0.8) : null;

      // Extract cropped region as separate base64 to boost accuracy
      const cropCanvas = document.createElement('canvas');
      cropCanvas.width = cropSide;
      cropCanvas.height = cropSide;
      const cropCtx = cropCanvas.getContext('2d');
      cropCtx.drawImage(canvas, cropX, cropY, cropSide, cropSide, 0, 0, cropSide, cropSide);
      const croppedImageDataUrl = cropCanvas.toDataURL('image/jpeg', 0.8);
      
      // Show transient outline at the clicked crop region (in container coords)
      try {
        const outlineSize = cropSide * scale;
        const displayCenterX = offsetX + centerX * scale;
        const displayCenterY = offsetY + centerY * scale;
        setClickPosition({
          x: Math.round(displayCenterX - outlineSize / 2),
          y: Math.round(displayCenterY - outlineSize / 2),
          size: Math.round(outlineSize)
        });
        setShowOutline(true);
        if (outlineTimeoutRef.current) {
          try { clearTimeout(outlineTimeoutRef.current); } catch (_) {}
        }
        outlineTimeoutRef.current = setTimeout(() => setShowOutline(false), 1000);
      } catch (_) {}

      // Multi-camera analysis: send both narrow and wide images if available
      let result;
      if (wideImageDataUrl) {
        // Send both images for better context
        result = await analyzeImage(
          fullImageDataUrl,
          Math.round(centerX),
          Math.round(centerY),
          cropX + Math.floor(cropSide / 2),
          cropY + Math.floor(cropSide / 2),
          cropSide,
          wideImageDataUrl // Pass wide image as additional context
        );
      } else {
        // Fallback to single camera AB test
        const abVariant = Math.random() < 0.5 ? 'coords' : 'crop';
        if (abVariant === 'coords') {
          result = await analyzeImage(
            fullImageDataUrl,
            Math.round(centerX),
            Math.round(centerY),
            cropX + Math.floor(cropSide / 2),
            cropY + Math.floor(cropSide / 2),
            cropSide
          );
        } else {
          // For cropped case, provide center relative to the crop
          const mid = Math.floor(cropSide / 2);
          result = await analyzeImage(
            croppedImageDataUrl,
            mid,
            mid,
            mid,
            mid,
            cropSide
          );
        }
      }
      if (onAnalysisComplete) {
        onAnalysisComplete({ ...result, hasMultipleViews: !!wideImageDataUrl });
      }
      } catch (error) {
        console.error('Camera structuralization failed:', error);
    } finally {
      setIsProcessing(false);
    }
  };

  const switchCamera = async () => {
    if (stream) {
      stream.getTracks().forEach(track => track.stop());
    }
    if (wideStream) {
      wideStream.getTracks().forEach(track => track.stop());
    }
    
    const currentFacingMode = stream?.getVideoTracks()[0]?.getSettings()?.facingMode || 'environment';
    const newFacingMode = currentFacingMode === 'environment' ? 'user' : 'environment';
    
    try {
      const mediaStream = await navigator.mediaDevices.getUserMedia({ 
        video: { 
          facingMode: newFacingMode,
          width: { ideal: 1920 },
          height: { ideal: 1080 },
          aspectRatio: { ideal: 16/9 }
        },
        audio: false 
      });
      
      if (videoRef.current) {
        videoRef.current.srcObject = mediaStream;
      }
      
      setStream(mediaStream);
      
      // Try to get wide camera for new facing mode
      try {
        const wideMediaStream = await navigator.mediaDevices.getUserMedia({
          video: {
            facingMode: newFacingMode,
            width: { ideal: 1920 },
            height: { ideal: 1080 },
            aspectRatio: { ideal: 16/9 },
            advanced: [{ zoom: { min: 0.5, max: 0.7 } }]
          },
          audio: false
        });
        
        if (wideVideoRef.current) {
          wideVideoRef.current.srcObject = wideMediaStream;
        }
        
        setWideStream(wideMediaStream);
      } catch (wideError) {
        console.log('Wide camera not available for', newFacingMode);
        setWideStream(null);
        setHasMultipleCameras(false);
      }
    } catch (error) {
      console.error('Camera switch error:', error);
    }
  };

  return (
    <div className="camera-container">
      <div
        className={`camera-box ${hasMultipleCameras ? 'multi-camera' : ''}`}
        onClick={handleCameraClick}
      >
        {/* Main camera preview */}
        <video
          ref={videoRef}
          autoPlay
          playsInline
          muted
          className="camera-video primary"
        />
        
        {/* Wide camera preview (if available) */}
        {hasMultipleCameras && (
          <video
            ref={wideVideoRef}
            autoPlay
            playsInline
            muted
            className="camera-video secondary"
          />
        )}
        
        {showOutline && clickPosition && (
          <div
            className="click-outline"
            style={{
              left: clickPosition.x,
              top: clickPosition.y,
              width: clickPosition.size,
              height: clickPosition.size
            }}
          />
        )}
        
        {hasMultipleCameras && (
          <div className="multi-camera-indicator">
            <span>📷 Multiple views</span>
          </div>
        )}
      </div>
      
      {showSwitchCamera && (
        <button onClick={switchCamera} className="switch-camera-btn">
          Switch Camera
        </button>
      )}
    </div>
  );
};

export default CameraSection;