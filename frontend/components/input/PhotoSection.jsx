import React, { useState } from 'react';
import { usePayment } from '../ui/PaymentContext';
import { useApi } from '../../hooks/useApi';

const PhotoSection = ({ isProcessing, setIsProcessing, setCurrentAnalysisType, onAnalysisComplete }) => {
  const [photoUrl, setPhotoUrl] = useState('');
  const { checkPaymentRequired } = usePayment();
  const { analyzeImage } = useApi();

  const handleFileUpload = async (e) => {
    const file = e.target.files?.[0];
    if (!file || isProcessing) return;

    if (checkPaymentRequired()) {
      return;
    }

    setIsProcessing(true);
    setCurrentAnalysisType('photo');

    try {
      // Convert file to base64
      const reader = new FileReader();
      reader.onload = async (event) => {
        const imageData = event.target.result;
        
        try {
          const result = await analyzeImage(imageData, file.name);
          if (onAnalysisComplete) {
            onAnalysisComplete(result);
          }
        } catch (error) {
          console.error('Photo analysis failed:', error);
        } finally {
          setIsProcessing(false);
        }
      };
      reader.readAsDataURL(file);
    } catch (error) {
      console.error('File reading failed:', error);
      setIsProcessing(false);
    }
  };

  const handleUrlAnalysis = async () => {
    if (!photoUrl.trim() || isProcessing) return;

    if (checkPaymentRequired()) {
      return;
    }

    setIsProcessing(true);
    setCurrentAnalysisType('url');

    try {
      // For URL, we send the URL itself to the backend
      const result = await analyzeImage(photoUrl, 'URL image');
      
      if (onAnalysisComplete) {
        onAnalysisComplete(result);
      }
      
      setPhotoUrl('');
    } catch (error) {
      console.error('URL analysis failed:', error);
    } finally {
      setIsProcessing(false);
    }
  };

  return (
    <div className="photo-options">
      <div className="upload-option">
        <input 
          type="file" 
          id="photo-upload" 
          accept="image/*"
          onChange={handleFileUpload}
          disabled={isProcessing}
        />
        <label htmlFor="photo-upload" className="upload-label">
          <svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
            <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"/>
            <polyline points="7,10 12,15 17,10"/>
            <line x1="12" y1="15" x2="12" y2="3"/>
          </svg>
          <span>Upload Photo</span>
        </label>
        
        <div className="url-input-container">
          <input 
            type="url" 
            id="photo-url" 
            placeholder="Paste image URL..."
            value={photoUrl}
            onChange={(e) => setPhotoUrl(e.target.value)}
            disabled={isProcessing}
          />
          <button 
            type="button" 
            className="url-button"
            onClick={handleUrlAnalysis}
            disabled={isProcessing || !photoUrl.trim()}
          >
            <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
              <circle cx="11" cy="11" r="8"/>
              <path d="M21 21l-4.35-4.35"/>
            </svg>
            <span>Analyze</span>
          </button>
        </div>
      </div>
    </div>
  );
};

export default PhotoSection;