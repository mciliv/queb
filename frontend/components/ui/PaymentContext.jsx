import React, { createContext, useContext, useState, useEffect } from 'react';
import logger from '../../core/logger.js';

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
  const [hasLoggedPaymentStatus, setHasLoggedPaymentStatus] = useState(false);

  useEffect(() => {
    if (!config.enabled) {
      if (!hasLoggedPaymentStatus) {
        logger.info('💳 Payment functionality disabled');
        setHasLoggedPaymentStatus(true);
      }
      return;
    }

    if (config.devMode) {
      if (!hasLoggedPaymentStatus) {
        logger.info('🔧 Auto-enabling developer mode for localhost');
        setHasLoggedPaymentStatus(true);
      }
      setIsPaymentSetup(true);
    }
  }, [config, hasLoggedPaymentStatus]);

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