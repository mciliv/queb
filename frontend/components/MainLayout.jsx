import React, { useState, useEffect, useRef } from 'react';
import TextInput from './TextInput';
import ModeSelector from './ModeSelector';
import CameraSection from './CameraSection';
import PhotoSection from './PhotoSection';
import Results from './Results';
import PaymentSection from './PaymentSection';
import { usePayment } from './PaymentContext';
import { useApi } from '../hooks/useApi';

const MainLayout = ({ 
  isProcessing,
  setIsProcessing,
  viewers,
  setViewers,
  currentAnalysisType,
  setCurrentAnalysisType,
  lastAnalysis,
  setLastAnalysis
}) => {
  const [objectInput, setObjectInput] = useState('');
  const [cameraMode, setCameraMode] = useState(false);
  const [photoMode, setPhotoMode] = useState(false);
  const { checkPaymentRequired } = usePayment();
  const { analyzeText } = useApi();

  // Handle keyboard shortcut
  useEffect(() => {
    const handleKeyDown = (event) => {
      if ((event.metaKey || event.ctrlKey) && event.key === 'k') {
        event.preventDefault();
        document.getElementById('object-input')?.focus();
      }
    };

    document.addEventListener('keydown', handleKeyDown);
    return () => document.removeEventListener('keydown', handleKeyDown);
  }, []);

  const handleTextAnalysis = async (value) => {
    if (isProcessing || !value.trim()) return;

    if (checkPaymentRequired()) {
      // Show payment modal
      return;
    }

    setIsProcessing(true);
    setCurrentAnalysisType('text');
    
    try {
      const result = await analyzeText(value);
      
      if (result.molecules && result.molecules.length > 0) {
        const newViewers = result.molecules.map(mol => ({
          name: mol.name || value,
          sdfData: mol.sdf_data,
          smiles: mol.smiles
        }));
        setViewers(prev => [...prev, ...newViewers]);
      }
      
      setLastAnalysis(result);
      setObjectInput('');
    } catch (error) {
      console.error('Analysis failed:', error);
    } finally {
      setIsProcessing(false);
    }
  };

  const handleAnalysisComplete = (result) => {
    if (result.molecules && result.molecules.length > 0) {
      const newViewers = result.molecules.map(mol => ({
        name: mol.name || 'Captured object',
        sdfData: mol.sdf_data,
        smiles: mol.smiles
      }));
      setViewers(prev => [...prev, ...newViewers]);
    }
    
    setLastAnalysis(result);
  };

  return (
    <div className="app-container">
      <div className="main-app-interface">
        <div className="main-content-layout">
          {/* Left side: Analysis section */}
          <div className="analysis-section">
            <TextInput 
              value={objectInput}
              onChange={setObjectInput}
              onSubmit={handleTextAnalysis}
              isProcessing={isProcessing}
            />

            <ModeSelector
              cameraMode={cameraMode}
              setCameraMode={setCameraMode}
              photoMode={photoMode}
              setPhotoMode={setPhotoMode}
            />

            {cameraMode && (
              <CameraSection
                isProcessing={isProcessing}
                setIsProcessing={setIsProcessing}
                setCurrentAnalysisType={setCurrentAnalysisType}
                onAnalysisComplete={handleAnalysisComplete}
              />
            )}

            {photoMode && (
              <PhotoSection
                isProcessing={isProcessing}
                setIsProcessing={setIsProcessing}
                setCurrentAnalysisType={setCurrentAnalysisType}
                onAnalysisComplete={handleAnalysisComplete}
              />
            )}

            <Results viewers={viewers} />
          </div>

          {/* Right side: Payment section */}
          <PaymentSection />
        </div>
      </div>
    </div>
  );
};

export default MainLayout;