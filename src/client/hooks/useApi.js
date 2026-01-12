/**
 * useApi Hook - Centralized API client for all backend communication
 * 
 * This custom React hook provides a clean interface for making API calls
 * with built-in features like caching, retry logic, and error handling.
 * It abstracts away the complexity of network requests and provides
 * consistent behavior across the application.
 * 
 * Key features:
 * - Request caching to avoid redundant API calls
 * - Automatic retry with exponential backoff
 * - Request deduplication (prevents duplicate concurrent requests)
 * - Loading and error state management
 * - Timeout handling
 * - Type-safe API methods for each endpoint
 * 
 * Usage:
 * const api = useApi();
 * const result = await api.structuralizeText('coffee', 'database');
 */

import { useState, useCallback, useRef } from 'react';

// Use same-origin for all API calls to avoid port/protocol mismatches
// This ensures the API works correctly in both development and production
const getApiBase = () => '';

const API_BASE = getApiBase();

// Configuration constants
const DEFAULT_RETRIES = 2;         // Retry failed requests twice
const RETRY_DELAY = 1000;          // 1 second initial delay between retries
const CACHE_DURATION = 300000;     // 5 minutes default cache
const MAX_CACHE_SIZE = 100;        // Prevent memory leaks

export const useApi = () => {
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null); // Full error object with diagnostics
  const requestCacheRef = useRef(new Map());
  const pendingRequestsRef = useRef(new Map());

  /**
   * Core API call function with caching, retry, and error handling
   * 
   * @param {string} endpoint - API endpoint path (e.g., '/api/structuralize')
   * @param {Object} options - Fetch options plus custom settings
   * @param {number} retryCount - Current retry attempt (used internally)
   * @returns {Promise<any>} - API response data
   * @throws {Error} - Network errors or API errors
   */
  const apiCall = useCallback(async (endpoint, options = {}, retryCount = 0) => {
    const maxRetries = options.maxRetries || DEFAULT_RETRIES;
    const cacheKey = `${endpoint}:${JSON.stringify(options.body || {})}`;
    const enableCaching = options.enableCaching !== false; // Default to true

    // Check cache first to avoid unnecessary API calls
    // Cache is used for GET requests or when explicitly enabled for POST
    if (enableCaching && (options.method !== 'POST' || options.cachePost)) {
      const cached = requestCacheRef.current.get(cacheKey);
      if (cached && Date.now() - cached.timestamp < (options.cacheDuration || CACHE_DURATION)) {
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

      // #region agent log
      fetch('http://127.0.0.1:7243/ingest/f1225f0b-6c5b-477f-bc5d-1e74641debf9',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'useApi.js:74',message:'fetch starting',data:{url:url,method:options.method||'GET',hasBody:!!options.body},timestamp:Date.now(),sessionId:'debug-load-failed',runId:'pre-fix',hypothesisId:'C'})}).catch(()=>{});
      // #endregion

      const response = await fetch(url, {
        headers: {
          'Content-Type': 'application/json',
          ...options.headers,
        },
        ...options,
      });

      // #region agent log
      fetch('http://127.0.0.1:7243/ingest/f1225f0b-6c5b-477f-bc5d-1e74641debf9',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'useApi.js:86',message:'fetch completed',data:{status:response.status,statusText:response.statusText,ok:response.ok},timestamp:Date.now(),sessionId:'debug-load-failed',runId:'pre-fix',hypothesisId:'C'})}).catch(()=>{});
      // #endregion

      if (!response.ok) {
        const errorData = await response.json().catch(() => ({}));
        const errorObj = {
          message: errorData.error || `API call failed: ${response.status} ${response.statusText}`,
          code: errorData.code,
          recoverable: errorData.recoverable,
          recovery: errorData.recovery,
          timestamp: errorData.timestamp,
          stack: errorData.stack,
          context: errorData.context,
          originalMessage: errorData.originalMessage
        };
        const err = new Error(errorObj.message);
        err.details = errorObj;
        throw err;
      }

      const data = await response.json();
      
      // Cache successful responses for performance
      if (enableCaching) {
        requestCacheRef.current.set(cacheKey, {
          data,
          timestamp: Date.now()
        });
        
        // Implement LRU cache eviction to prevent memory leaks
        // Remove oldest entry when cache exceeds size limit
        if (requestCacheRef.current.size > MAX_CACHE_SIZE) {
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
        errorMessage = `Request aborted (endpoint: ${endpoint}, attempt: ${retryCount + 1}/${maxRetries + 1})`;
        shouldRetry = retryCount < maxRetries;
      } else if (err.message.includes('Failed to fetch') || err.message.includes('NetworkError')) {
        errorMessage = `Network error: ${err.message} (endpoint: ${endpoint})`;
        shouldRetry = retryCount < maxRetries;
      } else if (err.message.includes('500') || err.message.includes('502') || err.message.includes('503')) {
        errorMessage = `Server error ${err.message.match(/\d{3}/)?.[0] || '5xx'}: ${endpoint}`;
        shouldRetry = retryCount < maxRetries;
      } else if (err.message.includes('429')) {
        errorMessage = `Rate limit exceeded (endpoint: ${endpoint})`;
        shouldRetry = retryCount < maxRetries;
      }

      // Retry logic
      if (shouldRetry) {
        await new Promise(resolve => setTimeout(resolve, RETRY_DELAY * (retryCount + 1)));
        return apiCall(endpoint, options, retryCount + 1);
      }

      // Report image structuralization failures to backend logger (fire-and-forget). Why?
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
            message: `Image analysis failed: ${errorMessage} (endpoint: ${endpoint})`,
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

      const errorObj = {
        message: errorMessage,
        recoverable: shouldRetry,
        timestamp: new Date().toISOString(),
        endpoint: endpoint,  // Include endpoint for diagnostics
        function: endpoint.includes('structuralize') ? 'analyzeText' : 
                  endpoint.includes('generate-sdfs') ? 'generateSDFs' :
                  endpoint.includes('structuralize-image') ? 'analyzeImage' : 'unknown'
      };
      setError(errorObj);
      const error = new Error(errorMessage);
      error.details = errorObj;
      throw error;
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

  /**
   * Chemical compounds from object specified by text
   *
   * Takes any text description (e.g., "coffee", "aspirin", "plastic bottle")
   * and returns the chemical compounds present in that object.
   *
   * @param {string} text - Object/material description to analyze
   * @param {string} lookupMode - Analysis mode: 'GPT-5' (AI)
   * @returns {Promise<{object: string, chemicals: Array}>} Analysis results
   *
   * @example
   * const result = await structuralizeText('coffee', 'database');
   * // Returns: { object: 'coffee', chemicals: [{name: 'Caffeine', sdfPath: '...'}] }
   */
  const structuralizeText = useCallback(async (text, lookupMode = 'GPT-5') => {
    // Always trim to remove whitespace
    const trimmedText = text.trim();

    // #region agent log
    fetch('http://127.0.0.1:7243/ingest/f1225f0b-6c5b-477f-bc5d-1e74641debf9',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'useApi.js:222',message:'structuralizeText called',data:{text:trimmedText,lookupMode:lookupMode,textLength:trimmedText.length},timestamp:Date.now(),sessionId:'debug-load-failed',runId:'pre-fix',hypothesisId:'A,C'})}).catch(()=>{});
    // #endregion

    try {
      const result = await apiCall('/api/structuralize', {
        method: 'POST',
        body: JSON.stringify({ text: trimmedText, lookupMode }),
        maxRetries: 2,
        cachePost: true,              // Cache results to avoid re-analyzing same text
        cacheDuration: 600000,        // 10 minutes - chemical data doesn't change often
      });

      // #region agent log
      fetch('http://127.0.0.1:7243/ingest/f1225f0b-6c5b-477f-bc5d-1e74641debf9',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'useApi.js:235',message:'structuralizeText success',data:{hasResult:!!result,resultKeys:Object.keys(result||{}).join(',')},timestamp:Date.now(),sessionId:'debug-load-failed',runId:'pre-fix',hypothesisId:'A,C'})}).catch(()=>{});
      // #endregion

      return result;
    } catch (error) {
      // #region agent log
      fetch('http://127.0.0.1:7243/ingest/f1225f0b-6c5b-477f-bc5d-1e74641debf9',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'useApi.js:242',message:'structuralizeText error',data:{error:error.message,hasDetails:!!error.details},timestamp:Date.now(),sessionId:'debug-load-failed',runId:'pre-fix',hypothesisId:'A,C'})}).catch(()=>{});
      // #endregion
      throw error;
    }
  }, [apiCall]);

  /**
   * Analyze image to identify object and its chemical compounds
   * 
   * Uses AI vision to identify what's in the image (especially at click coordinates)
   * then analyzes the chemical composition of the identified object.
   * 
   * @param {string} imageData - Base64 encoded image or data URL
   * @param {number} x - X coordinate where user clicked
   * @param {number} y - Y coordinate where user clicked  
   * @param {number} cropMiddleX - Center X of crop region
   * @param {number} cropMiddleY - Center Y of crop region
   * @param {number} cropSize - Size of square crop region
   * @returns {Promise<{object: string, chemicals: Array}>} Analysis results
   * 
   * @example
   * const result = await structuralizeImage(base64Image, 150, 200);
   * // AI identifies object at coordinates and returns its chemicals
   */
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

    // Clamp coordinates to reasonable bounds to prevent errors
    const clamp = (n, min, max) => Math.min(max, Math.max(min, n));
    const payload = { imageBase64: base64 };
    
    // Add coordinates if provided (helps AI focus on specific area)
    if (typeof x === 'number') payload.x = clamp(Math.round(x), 0, 1000);
    if (typeof y === 'number') payload.y = clamp(Math.round(y), 0, 1000);
    if (typeof cropMiddleX === 'number') payload.cropMiddleX = clamp(Math.round(cropMiddleX), 0, 1000);
    if (typeof cropMiddleY === 'number') payload.cropMiddleY = clamp(Math.round(cropMiddleY), 0, 1000);
    if (typeof cropSize === 'number') payload.cropSize = clamp(Math.round(cropSize), 10, 500);


    return apiCall('/api/structuralize-image', {
      method: 'POST',
      body: JSON.stringify(payload),
      maxRetries: 1,
    });
  }, [apiCall]);

  const generateSDFs = useCallback(async (smilesArray, overwrite = false) => {
    if (!smilesArray || !Array.isArray(smilesArray) || smilesArray.length === 0) {
      throw new Error('SMILES array is required');
    }

    return apiCall('/api/generate-sdfs', {
      method: 'POST',
      body: JSON.stringify({ 
        smiles: smilesArray,
        overwrite: overwrite
      }),
      maxRetries: 1,
      cachePost: !overwrite, // Cache if not overwriting
      cacheDuration: 1800000, // 30 minutes for SDF generation
    });
  }, [apiCall]);


  const nameToSdf = useCallback(async (name, overwrite = false) => {
    if (typeof name !== 'string' || name.trim().length === 0) {
      throw new Error('Valid name is required');
    }
    return apiCall('/api/name-to-sdf', {
      method: 'POST',
      body: JSON.stringify({ name, overwrite }),
      maxRetries: 1,
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
    clearError,
    clearCache,
    testConnection,
  };
};