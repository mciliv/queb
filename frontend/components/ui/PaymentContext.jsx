import React, { createContext, useContext, useState, useEffect } from 'react';

const PaymentContext = createContext();

export const usePayment = () => {
  const context = useContext(PaymentContext);
  if (!context) {
    throw new Error('usePayment must be used within PaymentProvider');
  }
  return context;
};

export const PaymentProvider = ({ children, config }) => {
  const [isPaymentSetup, setIsPaymentSetup] = useState(!config.enabled);
  const [paymentMethods, setPaymentMethods] = useState([]);
  const [selectedMethod, setSelectedMethod] = useState(null);
  const [isLoading, setIsLoading] = useState(false);

  useEffect(() => {
    if (!config.enabled) {
      console.log('💳 Payment functionality disabled');
      return;
    }

    if (config.devMode) {
      console.log('🔧 Auto-enabling developer mode for dev.queb.space');
      setIsPaymentSetup(true);
    }
  }, [config]);

  const checkPaymentRequired = () => {
    if (config.enabled && !isPaymentSetup && !config.devMode) {
      return true;
    }
    return false;
  };

  const value = {
    isPaymentSetup,
    setIsPaymentSetup,
    paymentMethods,
    setPaymentMethods,
    selectedMethod,
    setSelectedMethod,
    isLoading,
    setIsLoading,
    checkPaymentRequired,
    config
  };

  return (
    <PaymentContext.Provider value={value}>
      {children}
    </PaymentContext.Provider>
  );
};