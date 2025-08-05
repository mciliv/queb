import { useState, useCallback } from 'react';

const API_BASE = window.location.hostname === 'localhost' ? 'https://localhost:3002' : '';

export const useApi = () => {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const apiCall = useCallback(async (endpoint, options = {}) => {
    setLoading(true);
    setError(null);

    try {
      const url = `${API_BASE}${endpoint}`;
      const response = await fetch(url, {
        headers: {
          'Content-Type': 'application/json',
          ...options.headers,
        },
        ...options,
      });

      if (!response.ok) {
        throw new Error(`API call failed: ${response.statusText}`);
      }

      const data = await response.json();
      return data;
    } catch (err) {
      setError(err.message);
      throw err;
    } finally {
      setLoading(false);
    }
  }, []);

  const analyzeText = useCallback(async (text) => {
    return apiCall('/analyze-text', {
      method: 'POST',
      body: JSON.stringify({ object: text }),
    });
  }, [apiCall]);

  const analyzeImage = useCallback(async (imageData, objectName) => {
    return apiCall('/analyze-image', {
      method: 'POST',
      body: JSON.stringify({ 
        image: imageData,
        object_name: objectName 
      }),
    });
  }, [apiCall]);

  return {
    loading,
    error,
    apiCall,
    analyzeText,
    analyzeImage,
  };
};