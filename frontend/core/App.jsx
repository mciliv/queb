import React, { useState, useEffect, useRef } from 'react';
import { PaymentProvider } from '../components/ui/PaymentContext';
import MainLayout from '../components/ui/MainLayout';
import '../assets/style.css';

// Payment toggle configuration
const PAYMENT_CONFIG = {
  enabled: false, // Set to true to enable payment functionality
  devMode: window.location.hostname === 'localhost' || window.location.hostname === '127.0.0.1'
};

function App() {
  const [isPaymentSetup, setIsPaymentSetup] = useState(!PAYMENT_CONFIG.enabled);
  const [isProcessing, setIsProcessing] = useState(false);
  const [viewers, setViewers] = useState([]);
  const [currentAnalysisType, setCurrentAnalysisType] = useState(null);
  const [lastAnalysis, setLastAnalysis] = useState(null);

  useEffect(() => {
    // Auto-enable dev mode for localhost
    if (PAYMENT_CONFIG.devMode) {
      console.log('🔧 Auto-enabling developer mode for localhost');
      setIsPaymentSetup(true);
    }
    
    console.log('✅ Molecular analysis app initialized');
  }, []);

  return (
    <PaymentProvider config={PAYMENT_CONFIG}>
      <MainLayout />
    </PaymentProvider>
  );
}

export default App;