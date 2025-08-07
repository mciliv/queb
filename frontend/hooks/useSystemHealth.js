import { useState, useEffect, useCallback } from 'react';

export const useSystemHealth = () => {
  const [healthResults, setHealthResults] = useState([]);
  const [errors, setErrors] = useState([]);
  const [isChecking, setIsChecking] = useState(false);

  const logResult = useCallback((test, status, message, data = null) => {
    const result = {
      test,
      status,
      message,
      timestamp: new Date().toISOString(),
      data
    };
    
    setHealthResults(prev => [...prev, result]);
    console.log(`[${status}] ${test}: ${message}`);
    
    if (data) {
      console.log('Data:', data);
    }
    
    return result;
  }, []);

  const checkDOMElements = useCallback(() => {
    const elements = [
      'video-feed',
      'object-input', 
      'payment-modal',
      'account-status'
    ];

    elements.forEach(id => {
      const element = document.getElementById(id);
      if (element) {
        logResult('DOM', 'PASS', `Element ${id} found`);
      } else {
        logResult('DOM', 'WARN', `Element ${id} not found`);
      }
    });
  }, [logResult]);

  const checkAPIEndpoints = useCallback(async () => {
    const endpoints = [
      { url: '/api/health', name: 'Health Check' },
      { url: '/api/analyze', name: 'Analysis API' },
      { url: '/api/user/check', name: 'User API' }
    ];

    for (const endpoint of endpoints) {
      try {
        const response = await fetch(endpoint.url);
        if (response.ok) {
          logResult('API', 'PASS', `${endpoint.name} responding`);
        } else {
          logResult('API', 'FAIL', `${endpoint.name} error: ${response.status}`);
        }
      } catch (error) {
        logResult('API', 'FAIL', `${endpoint.name} unreachable: ${error.message}`);
      }
    }
  }, [logResult]);

  const checkCameraSystem = useCallback(async () => {
    if (!navigator.mediaDevices || !navigator.mediaDevices.getUserMedia) {
      logResult('Camera', 'FAIL', 'getUserMedia not supported');
      return;
    }

    try {
      const stream = await navigator.mediaDevices.getUserMedia({ video: true });
      logResult('Camera', 'PASS', 'Camera access granted');
      stream.getTracks().forEach(track => track.stop());
    } catch (error) {
      logResult('Camera', 'WARN', `Camera access denied: ${error.message}`);
    }
  }, [logResult]);

  const checkPaymentSystem = useCallback(() => {
    if (typeof window.Stripe !== 'undefined') {
      logResult('Payment', 'PASS', 'Stripe loaded');
    } else {
      logResult('Payment', 'WARN', 'Stripe not loaded');
    }
  }, [logResult]);

  const checkPermissions = useCallback(async () => {
    if ('permissions' in navigator) {
      try {
        const cameraPermission = await navigator.permissions.query({ name: 'camera' });
        logResult('Permissions', 'INFO', `Camera permission: ${cameraPermission.state}`);
      } catch (error) {
        logResult('Permissions', 'WARN', `Cannot check camera permission: ${error.message}`);
      }
    }
  }, [logResult]);

  const runHealthCheck = useCallback(async () => {
    setIsChecking(true);
    setHealthResults([]);
    setErrors([]);
    
    console.log('🔍 Starting comprehensive health check...');
    
    checkDOMElements();
    checkPaymentSystem();
    await checkCameraSystem();
    await checkAPIEndpoints();
    await checkPermissions();
    
    setIsChecking(false);
    console.log('✅ Health check complete');
  }, [checkDOMElements, checkPaymentSystem, checkCameraSystem, checkAPIEndpoints, checkPermissions]);

  const generateReport = useCallback(() => {
    const passed = healthResults.filter(r => r.status === 'PASS').length;
    const failed = healthResults.filter(r => r.status === 'FAIL').length;
    const warnings = healthResults.filter(r => r.status === 'WARN').length;
    
    return {
      total: healthResults.length,
      passed,
      failed,
      warnings,
      score: healthResults.length > 0 ? (passed / healthResults.length) * 100 : 0
    };
  }, [healthResults]);

  return {
    healthResults,
    errors,
    isChecking,
    runHealthCheck,
    generateReport
  };
};