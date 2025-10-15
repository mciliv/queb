import { useState, useCallback, useRef } from 'react';

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
  const requestCacheRef = useRef(new Map());
  const pendingRequestsRef = useRef(new Map());

  const apiCall = useCallback(async (endpoint, options = {}, retryCount = 0) => {
    const maxRetries = options.maxRetries || DEFAULT_RETRIES;
    const timeout = options.timeout || DEFAULT_TIMEOUT;
    const cacheKey = `${endpoint}:${JSON.stringify(options.body || {})}`;
    const enableCaching = options.enableCaching !== false; // Default to true

    // Check cache first for GET requests or explicitly cached requests
    if (enableCaching && (options.method !== 'POST' || options.cachePost)) {
      const cached = requestCacheRef.current.get(cacheKey);
      if (cached && Date.now() - cached.timestamp < (options.cacheDuration || 300000)) { // 5 min default
        return cached.data;
      }
    }

    // Check for pending identical requests to avoid duplicates
    if (pendingRequestsRef.current.has(cacheKey)) {
      return pendingRequestsRef.current.get(cacheKey);
    }

    setLoading(true);
    setError(null);

    const requestPromise = (async () => {
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
      
      // Cache successful responses
      if (enableCaching) {
        requestCacheRef.current.set(cacheKey, {
          data,
          timestamp: Date.now()
        });
        // Limit cache size to prevent memory leaks
        if (requestCacheRef.current.size > 100) {
          const firstKey = requestCacheRef.current.keys().next().value;
          requestCacheRef.current.delete(firstKey);
        }
      }
      
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

      // Report image structuralization failures to backend logger (fire-and-forget)
      try {
        const isImageStructuralization = typeof endpoint === 'string' && (
          endpoint.includes('structures-from-image') ||
          endpoint.includes('image-molecules') ||
          endpoint.includes('estimate-image')
        );
        if (isImageStructuralization) {
          const payload = {
            timestamp: new Date().toISOString(),
            type: 'error',
            source: 'frontend',
            message: `AI failed to analyze image: ${errorMessage}`,
            endpoint,
          };
          if (typeof navigator !== 'undefined' && navigator.sendBeacon) {
            const blob = new Blob([JSON.stringify(payload)], { type: 'application/json' });
            navigator.sendBeacon('/api/log-error', blob);
          } else {
            // Best-effort logging without blocking UI
            fetch('/api/log-error', {
              method: 'POST',
              headers: { 'Content-Type': 'application/json' },
              body: JSON.stringify(payload),
              keepalive: true,
            }).catch(() => {});
          }
        }
      } catch (_) {
        // ignore reporter failures
      }

      setError(errorMessage);
      throw new Error(errorMessage);
    } finally {
      setLoading(false);
      pendingRequestsRef.current.delete(cacheKey);
    }
    })();
    
    // Store pending request
    pendingRequestsRef.current.set(cacheKey, requestPromise);
    
    try {
      const result = await requestPromise;
      return result;
    } catch (error) {
      pendingRequestsRef.current.delete(cacheKey);
      throw error;
    }
  }, []);

  const structuralizeText = useCallback(async (text, lookupMode = 'database') => {
    if (!text || !text.trim()) {
      throw new Error('Text input is required');
    }

    return apiCall('/structuralize', {
      method: 'POST',
      body: JSON.stringify({ object: text, lookupMode }),
      maxRetries: 2,
      timeout: 30000,
      cachePost: true,
      cacheDuration: 600000, // 10 minutes for text analysis
    });
  }, [apiCall]);

  const structuralizeImage = useCallback(async (
    imageData,
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

    const clamp = (n, min, max) => Math.min(max, Math.max(min, n));
    const payload = { imageBase64: base64 };
    if (typeof x === 'number') payload.x = clamp(Math.round(x), 0, 1000);
    if (typeof y === 'number') payload.y = clamp(Math.round(y), 0, 1000);
    if (typeof cropMiddleX === 'number') payload.cropMiddleX = clamp(Math.round(cropMiddleX), 0, 1000);
    if (typeof cropMiddleY === 'number') payload.cropMiddleY = clamp(Math.round(cropMiddleY), 0, 1000);
    if (typeof cropSize === 'number') payload.cropSize = clamp(Math.round(cropSize), 10, 500);

    try {
      const previewType = typeof imageData === 'string' && imageData.startsWith('data:')
        ? imageData.slice(5, imageData.indexOf(';'))
        : 'base64';
      // Structured, safe debug info (no base64 dump)
      console.log('AnalyzeImage payload summary:', {
        type: previewType,
        hasCoords: typeof x === 'number' && typeof y === 'number',
        coords: typeof x === 'number' && typeof y === 'number' ? { x: payload.x, y: payload.y } : null,
        hasCrop: typeof cropMiddleX === 'number' && typeof cropMiddleY === 'number' && typeof cropSize === 'number',
        crop: (typeof cropMiddleX === 'number' && typeof cropMiddleY === 'number' && typeof cropSize === 'number')
          ? { cx: payload.cropMiddleX, cy: payload.cropMiddleY, size: payload.cropSize }
          : null
      });
    } catch (_) {}

    return apiCall('/structuralize', {
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
      timeout: 30000,
      cachePost: !overwrite, // Cache if not overwriting
      cacheDuration: 1800000, // 30 minutes for SDF generation
    });
  }, [apiCall]);


  const nameToSdf = useCallback(async (name, overwrite = false) => {
    if (typeof name !== 'string' || name.trim().length === 0) {
      throw new Error('Valid name is required');
    }
    return apiCall('/name-to-sdf', {
      method: 'POST',
      body: JSON.stringify({ name, overwrite }),
      maxRetries: 1,
      timeout: 20000,
      cachePost: !overwrite,
      cacheDuration: 1800000,
    });
  }, [apiCall]);

  const clearError = useCallback(() => {
    setError(null);
  }, []);

  const clearCache = useCallback(() => {
    requestCacheRef.current.clear();
    pendingRequestsRef.current.clear();
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
    analyzeText: structuralizeText,
    structuresFromText: structuralizeText,
    analyzeImage: structuralizeImage,
    generateSDFs,
    nameToSdf,
    // FoodDB methods for input methods to call
    listFoods: useCallback(async (limit = 25) => {
      const params = new URLSearchParams({ limit: String(limit) });
      return apiCall(`/api/fooddb/foods?${params.toString()}`, { method: 'GET' });
    }, [apiCall]),
    searchFoods: useCallback(async (query) => {
      const q = typeof query === 'string' ? query.trim() : '';
      const params = new URLSearchParams({ q });
      return apiCall(`/api/fooddb/search?${params.toString()}`, { method: 'GET' });
    }, [apiCall]),
    getFood: useCallback(async (id) => {
      return apiCall(`/api/fooddb/foods/${encodeURIComponent(String(id))}`, { method: 'GET' });
    }, [apiCall]),
    getFoodCompounds: useCallback(async (id, name = null) => {
      const params = name ? `?name=${encodeURIComponent(name)}` : '';
      return apiCall(`/api/fooddb/foods/${encodeURIComponent(String(id))}/compounds${params}`, { method: 'GET' });
    }, [apiCall]),
    clearError,
    clearCache,
    testConnection,
  };
};