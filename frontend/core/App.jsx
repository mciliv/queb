import React, { useState, useEffect, useRef } from 'react';
import { PaymentProvider } from '../components/PaymentContext';
import MainLayout from '../components/MainLayout';
import '../assets/style.css';

// Payment toggle configuration
const PAYMENT_CONFIG = {
  enabled: false, // Set to true to enable payment functionality
  devMode: window.location.hostname === 'dev.queb.space'
};

function App() {
  const [isPaymentSetup, setIsPaymentSetup] = useState(!PAYMENT_CONFIG.enabled);
  const [isProcessing, setIsProcessing] = useState(false);
  const [viewers, setViewers] = useState([]);
  const [currentAnalysisType, setCurrentAnalysisType] = useState(null);
  const [lastAnalysis, setLastAnalysis] = useState(null);

  useEffect(() => {
    // Auto-enable dev mode for dev.queb.space
    if (PAYMENT_CONFIG.devMode) {
      console.log('🔧 Auto-enabling developer mode for dev.queb.space');
      setIsPaymentSetup(true);
    }
    
    console.log('✅ Molecular analysis app initialized');
  }, []);

  return (
    <PaymentProvider config={PAYMENT_CONFIG}>
      <MainLayout 
        isProcessing={isProcessing}
        setIsProcessing={setIsProcessing}
        viewers={viewers}
        setViewers={setViewers}
        currentAnalysisType={currentAnalysisType}
        setCurrentAnalysisType={setCurrentAnalysisType}
        lastAnalysis={lastAnalysis}
        setLastAnalysis={setLastAnalysis}
        isPaymentSetup={isPaymentSetup}
        setIsPaymentSetup={setIsPaymentSetup}
      />
    </PaymentProvider>
  );
}

export default App;