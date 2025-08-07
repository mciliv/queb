import { useState, useEffect, useCallback, useRef } from 'react';

export const useUIContext = () => {
  const [snapshots, setSnapshots] = useState(new Map());
  const [changeHistory, setChangeHistory] = useState([]);
  const [enabled, setEnabled] = useState(false);
  const observerRef = useRef(null);

  const generateUISnapshot = useCallback(() => {
    return {
      timestamp: new Date().toISOString(),
      viewport: {
        width: window.innerWidth,
        height: window.innerHeight,
        devicePixelRatio: window.devicePixelRatio
      },
      layout: {
        sidebarVisible: !document.getElementById('payment-section')?.classList.contains('hidden'),
        accountButtonPosition: getElementPosition('.account-link'),
        textInputMargin: getComputedStyle(document.getElementById('object-input'))?.marginRight,
        safeAreaRight: getComputedStyle(document.documentElement).getPropertyValue('--safe-area-right')
      },
      activeElements: {
        focusedElement: document.activeElement?.id || null,
        visibleSections: getVisibleSections(),
        errorStates: getErrorStates()
      },
      cssVariables: getCSSVariables(),
      responsiveState: getResponsiveBreakpoint()
    };
  }, []);

  const getElementPosition = useCallback((selector) => {
    const element = document.querySelector(selector);
    if (!element) return null;
    const rect = element.getBoundingClientRect();
    return {
      top: rect.top,
      left: rect.left,
      width: rect.width,
      height: rect.height
    };
  }, []);

  const getVisibleSections = useCallback(() => {
    const sections = ['main-content', 'payment-section', 'account-modal'];
    return sections.filter(id => {
      const element = document.getElementById(id);
      return element && !element.classList.contains('hidden');
    });
  }, []);

  const getErrorStates = useCallback(() => {
    const errorElements = document.querySelectorAll('.error, .alert, .warning');
    return Array.from(errorElements).map(el => ({
      id: el.id || 'unnamed',
      message: el.textContent?.trim(),
      visible: !el.classList.contains('hidden')
    }));
  }, []);

  const getCSSVariables = useCallback(() => {
    const style = getComputedStyle(document.documentElement);
    return {
      primaryColor: style.getPropertyValue('--primary-color'),
      backgroundColor: style.getPropertyValue('--background-color'),
      textColor: style.getPropertyValue('--text-color'),
      safeAreaRight: style.getPropertyValue('--safe-area-right')
    };
  }, []);

  const getResponsiveBreakpoint = useCallback(() => {
    const width = window.innerWidth;
    if (width < 768) return 'mobile';
    if (width < 1024) return 'tablet';
    return 'desktop';
  }, []);

  const trackChange = useCallback((type, details) => {
    const change = {
      type,
      details,
      timestamp: new Date().toISOString(),
      snapshot: generateUISnapshot()
    };
    
    setChangeHistory(prev => [...prev.slice(-49), change]); // Keep last 50 changes
    console.log(`🔄 UI Change tracked: ${type}`, details);
  }, [generateUISnapshot]);

  const setupMutationObserver = useCallback(() => {
    if (observerRef.current) {
      observerRef.current.disconnect();
    }

    observerRef.current = new MutationObserver((mutations) => {
      mutations.forEach((mutation) => {
        if (mutation.type === 'attributes' && mutation.attributeName === 'class') {
          trackChange('class-change', {
            element: mutation.target.id || mutation.target.className,
            newValue: mutation.target.getAttribute('class')
          });
        }
        
        if (mutation.type === 'childList' && mutation.addedNodes.length > 0) {
          trackChange('dom-addition', {
            parentElement: mutation.target.id || mutation.target.className,
            addedNodes: mutation.addedNodes.length
          });
        }
      });
    });

    observerRef.current.observe(document.body, {
      attributes: true,
      childList: true,
      subtree: true,
      attributeFilter: ['class', 'style', 'hidden']
    });
  }, [trackChange]);

  const enable = useCallback(() => {
    setEnabled(true);
    setupMutationObserver();
    
    // Track resize events
    const handleResize = () => trackChange('viewport-resize', {
      width: window.innerWidth,
      height: window.innerHeight
    });
    window.addEventListener('resize', handleResize);
    
    console.log('🤖 UI Context Automation enabled');
    
    return () => {
      window.removeEventListener('resize', handleResize);
      if (observerRef.current) {
        observerRef.current.disconnect();
      }
    };
  }, [setupMutationObserver, trackChange]);

  const disable = useCallback(() => {
    setEnabled(false);
    if (observerRef.current) {
      observerRef.current.disconnect();
    }
    console.log('🤖 UI Context Automation disabled');
  }, []);

  const captureSnapshot = useCallback((label = 'manual') => {
    const snapshot = generateUISnapshot();
    setSnapshots(prev => new Map(prev.set(label, snapshot)));
    return snapshot;
  }, [generateUISnapshot]);

  // Cleanup on unmount
  useEffect(() => {
    return () => {
      if (observerRef.current) {
        observerRef.current.disconnect();
      }
    };
  }, []);

  return {
    enabled,
    snapshots,
    changeHistory,
    enable,
    disable,
    captureSnapshot,
    generateUISnapshot,
    trackChange
  };
};