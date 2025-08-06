import { useState, useCallback } from 'react';

const API_BASE = window.location.hostname === 'localhost' ? 'http://localhost:8080' : '';

// Default timeout and retry settings
const DEFAULT_TIMEOUT = 30000; // 30 seconds
const DEFAULT_RETRIES = 2;
const RETRY_DELAY = 1000; // 1 second

export const useApi = () => {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const apiCall = useCallback(async (endpoint, options = {}, retryCount = 0) => {
    const maxRetries = options.maxRetries || DEFAULT_RETRIES;
    const timeout = options.timeout || DEFAULT_TIMEOUT;

    setLoading(true);
    setError(null);

    try {
      const url = `${API_BASE}${endpoint}`;
      
      // Create abort controller for timeout
      const controller = new AbortController();
      const timeoutId = setTimeout(() => controller.abort(), timeout);

      const response = await fetch(url, {
        headers: {
          'Content-Type': 'application/json',
          ...options.headers,
        },
        signal: controller.signal,
        ...options,
      });

      clearTimeout(timeoutId);

      if (!response.ok) {
        const errorData = await response.json().catch(() => ({}));
        throw new Error(errorData.error || `API call failed: ${response.status} ${response.statusText}`);
      }

      const data = await response.json();
      return data;
    } catch (err) {
      // Handle different types of errors
      let errorMessage = err.message;
      let shouldRetry = false;

      if (err.name === 'AbortError') {
        errorMessage = 'Request timed out. Please try again.';
        shouldRetry = retryCount < maxRetries;
      } else if (err.message.includes('Failed to fetch') || err.message.includes('NetworkError')) {
        errorMessage = 'Network connection failed. Please check your internet connection.';
        shouldRetry = retryCount < maxRetries;
      } else if (err.message.includes('500') || err.message.includes('502') || err.message.includes('503')) {
        errorMessage = 'Server error. Please try again in a moment.';
        shouldRetry = retryCount < maxRetries;
      } else if (err.message.includes('429')) {
        errorMessage = 'Rate limit exceeded. Please wait a moment before trying again.';
        shouldRetry = retryCount < maxRetries;
      }

      // Retry logic
      if (shouldRetry) {
        console.log(`API call failed, retrying... (${retryCount + 1}/${maxRetries + 1})`);
        await new Promise(resolve => setTimeout(resolve, RETRY_DELAY * (retryCount + 1)));
        return apiCall(endpoint, options, retryCount + 1);
      }

      setError(errorMessage);
      throw new Error(errorMessage);
    } finally {
      setLoading(false);
    }
  }, []);

  const analyzeText = useCallback(async (text) => {
    if (!text || !text.trim()) {
      throw new Error('Text input is required');
    }

    return apiCall('/analyze-text', {
      method: 'POST',
      body: JSON.stringify({ object: text }),
      maxRetries: 2,
      timeout: 45000, // 45 seconds for text analysis
    });
  }, [apiCall]);

  const analyzeImage = useCallback(async (imageData, objectName) => {
    if (!imageData) {
      throw new Error('Image data is required');
    }

    return apiCall('/analyze-image', {
      method: 'POST',
      body: JSON.stringify({ 
        image: imageData,
        object_name: objectName 
      }),
      maxRetries: 1, // Don't retry image analysis as much
      timeout: 60000, // 60 seconds for image analysis
    });
  }, [apiCall]);

  const generateSDFs = useCallback(async (smilesArray, overwrite = false) => {
    if (!smilesArray || !Array.isArray(smilesArray) || smilesArray.length === 0) {
      throw new Error('SMILES array is required');
    }

    return apiCall('/generate-sdfs', {
      method: 'POST',
      body: JSON.stringify({ 
        smiles: smilesArray,
        overwrite: overwrite
      }),
      maxRetries: 1,
      timeout: 30000, // 30 seconds for SDF generation
    });
  }, [apiCall]);

  const clearError = useCallback(() => {
    setError(null);
  }, []);

  return {
    loading,
    error,
    apiCall,
    analyzeText,
    analyzeImage,
    generateSDFs,
    clearError,
  };
};