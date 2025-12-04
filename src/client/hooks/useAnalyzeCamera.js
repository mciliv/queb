/**
 * useAnalyzeCamera - Live camera prediction for molecular discovery
 * 
 * This hook provides a focused interface for camera-based prediction.
 * Users can point their camera at objects and click to predict their molecular composition.
 * 
 * @example
 * const { 
 *   startCamera, 
 *   stopCamera, 
 *   analyzeClick,
 *   isAnalyzing 
 * } = useAnalyzeCamera();
 * 
 * // User clicks on coffee cup in camera view
 * const handleCameraClick = async (x, y) => {
 *   const result = await analyzeClick(x, y);
 *   // Shows molecules in the coffee
 * };
 */

import { useState, useCallback, useRef } from 'react';
import { useApi } from './useApi';  // Compose on top of existing hook

export function useAnalyzeCamera() {
  const [isAnalyzing, setIsAnalyzing] = useState(false);
  const [molecules, setMolecules] = useState([]);
  const [identifiedObject, setIdentifiedObject] = useState('');
  const [error, setError] = useState(null);
  const [cameraActive, setCameraActive] = useState(false);
  
  const videoRef = useRef(null);
  const streamRef = useRef(null);
  const { structuralizeImage } = useApi();  // Use existing functionality
  
  /**
   * Start camera stream
   * @param {HTMLVideoElement} videoElement - Video element to stream to
   * @returns {Promise<MediaStream>} Camera stream
   */
  const startCamera = useCallback(async (videoElement) => {
    try {
      const stream = await navigator.mediaDevices.getUserMedia({
        video: {
          facingMode: 'environment', // Prefer back camera
          width: { ideal: 1920 },
          height: { ideal: 1080 }
        },
        audio: false
      });
      
      videoRef.current = videoElement;
      streamRef.current = stream;
      
      if (videoElement) {
        videoElement.srcObject = stream;
      }
      
      setCameraActive(true);
      return stream;
      
    } catch (err) {
      const message = 'Camera access denied or unavailable';
      setError(message);
      throw new Error(message);
    }
  }, []);
  
  /**
   * Stop camera stream
   */
  const stopCamera = useCallback(() => {
    if (streamRef.current) {
      streamRef.current.getTracks().forEach(track => track.stop());
      streamRef.current = null;
    }
    
    if (videoRef.current) {
      videoRef.current.srcObject = null;
      videoRef.current = null;
    }
    
    setCameraActive(false);
  }, []);
  
  /**
   * Capture frame from video
   */
  const captureFrame = useCallback(() => {
    if (!videoRef.current) {
      throw new Error('Camera not active');
    }
    
    const video = videoRef.current;
    const canvas = document.createElement('canvas');
    canvas.width = video.videoWidth;
    canvas.height = video.videoHeight;
    
    const ctx = canvas.getContext('2d');
    ctx.drawImage(video, 0, 0);
    
    return canvas.toDataURL('image/jpeg');
  }, []);
  
  /**
   * Analyze what user clicked on in camera view
   * @param {number} x - Click X coordinate in video
   * @param {number} y - Click Y coordinate in video
   * @param {Object} options - Analysis options
   * @returns {Promise<{object: string, molecules: Array}>} Analysis results
   */
  const analyzeClick = useCallback(async (x, y, options = {}) => {
    if (!cameraActive || !videoRef.current) {
      throw new Error('Camera must be active to analyze');
    }
    
    setIsAnalyzing(true);
    setError(null);
    
    try {
      // Capture current frame
      const imageData = captureFrame();
      const base64 = imageData.split(',')[1];
      
      // Calculate crop region (25% of smallest dimension)
      const video = videoRef.current;
      const cropSize = Math.min(video.videoWidth, video.videoHeight) * 0.25;
      
      // Use existing API method
      const result = await structuralizeImage(
        imageData,  // full data URL
        Math.round(x),
        Math.round(y),
        Math.round(x),  // cropMiddleX
        Math.round(y),  // cropMiddleY
        Math.round(cropSize)
      );
      
      // Transform to user-friendly format
      const molecules = (result.chemicals || result.molecules || []).map(chem => ({
        name: chem.name || chem.chemicalName,
        structure3D: chem.sdfPath,
        status: chem.status || 'ok'
      }));
      
      setMolecules(molecules);
      setIdentifiedObject(result.object || result.result?.object || 'Unknown object');
      
      return {
        object: result.object,
        molecules,
        clickPoint: { x, y },
        cropRegion: {
          x: x - cropSize/2,
          y: y - cropSize/2,
          width: cropSize,
          height: cropSize
        }
      };
      
    } catch (err) {
      const message = err.message || 'Failed to analyze camera view';
      setError(message);
      throw new Error(message);
    } finally {
      setIsAnalyzing(false);
    }
  }, [cameraActive, captureFrame, post]);
  
  /**
   * Clear results
   */
  const clear = useCallback(() => {
    setMolecules([]);
    setIdentifiedObject('');
    setError(null);
  }, []);
  
  return {
    // Camera controls
    startCamera,
    stopCamera,
    
    // Analysis
    analyzeClick,
    captureFrame,
    clear,
    
    // State
    isAnalyzing,
    cameraActive,
    molecules,
    identifiedObject,
    error,
    
    // Helpers
    hasMolecules: molecules.length > 0,
    canAnalyze: cameraActive && !isAnalyzing
  };
}


