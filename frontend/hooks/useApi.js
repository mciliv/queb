import { useState, useCallback } from 'react';

// Use same-origin for all API calls to avoid port/protocol mismatches
const getApiBase = () => '';

const API_BASE = getApiBase();

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
        await new Promise(resolve => setTimeout(resolve, RETRY_DELAY * (retryCount + 1)));
        return apiCall(endpoint, options, retryCount + 1);
      }

      setError(errorMessage);
      throw new Error(errorMessage);
    } finally {
      setLoading(false);
    }
  }, []);

  const structuresFromText = useCallback(async (text) => {
    if (!text || !text.trim()) {
      throw new Error('Text input is required');
    }

    return apiCall('/structures-from-text', {
      method: 'POST',
      body: JSON.stringify({ object: text }),
      maxRetries: 2,
      timeout: 30000,
    });
  }, [apiCall]);

  const structuralizeImage = useCallback(async (
    imageData,
    _label,
    x,
    y,
    cropMiddleX,
    cropMiddleY,
    cropSize
  ) => {
    if (!imageData) {
      throw new Error('Image data is required');
    }

    // Strip data URL prefix if present; backend expects raw base64
    const base64 = typeof imageData === 'string' && imageData.startsWith('data:')
      ? imageData.split(',')[1]
      : imageData;

    const payload = { imageBase64: base64 };
    if (typeof x === 'number') payload.x = x;
    if (typeof y === 'number') payload.y = y;
    if (typeof cropMiddleX === 'number') payload.cropMiddleX = cropMiddleX;
    if (typeof cropMiddleY === 'number') payload.cropMiddleY = cropMiddleY;
    if (typeof cropSize === 'number') payload.cropSize = cropSize;

    return apiCall('/structuralize-image', {
      method: 'POST',
      body: JSON.stringify(payload),
      maxRetries: 1,
      timeout: 60000,
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

  const twoNamesToSdf = useCallback(async (names, overwrite = false) => {
    let list = names;
    if (typeof names === 'string') {
      list = names.split(',').map(s => s.trim()).filter(Boolean);
    }
    if (!Array.isArray(list) || list.length === 0) {
      throw new Error('Provide one or two names');
    }
    return apiCall('/two-names-to-sdf', {
      method: 'POST',
      body: JSON.stringify({ names: list, overwrite }),
      maxRetries: 1,
      timeout: 30000,
    });
  }, [apiCall]);

  const clearError = useCallback(() => {
    setError(null);
  }, []);

  const testConnection = useCallback(async () => {
    try {
      const url = `${API_BASE}/health`;
      const response = await fetch(url);
      return response.ok;
    } catch (_) {
      return false;
    }
  }, []);

  return {
    loading,
    error,
    apiCall,
    analyzeText: structuresFromText,
    structuresFromText,
    analyzeImage: structuralizeImage,
    generateSDFs,
    twoNamesToSdf,
    clearError,
    testConnection,
  };
};