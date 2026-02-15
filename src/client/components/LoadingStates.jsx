import React from 'react';
import { useAppController } from '../AppController';
import '../assets/style.css';

const LoadingStates = () => {
  const { isProcessing } = useAppController();

  if (!isProcessing) {
    return null;
  }

  return (
    <div className="loading-overlay">
      <div className="spinner"></div>
    </div>
  );
};

export default LoadingStates;