import React, { useState } from 'react';
import { useApi } from '../../hooks/useApi';

const LinkSection = ({ isProcessing, setIsProcessing, onAnalysisComplete }) => {
  const [linkInput, setLinkInput] = useState('');
  const [error, setError] = useState('');
  const { analyzeImage } = useApi();

  const handleLinkSubmit = async (e) => {
    e.preventDefault();
    
    if (!linkInput.trim()) {
      setError('Please enter an image URL');
      return;
    }

    // Basic URL validation
    try {
      new URL(linkInput);
    } catch {
      setError('Please enter a valid URL');
      return;
    }

    setIsProcessing(true);
    setError('');

    try {
      // Pass URL as image data - backend will handle fetching from URL
      const result = await analyzeImage(linkInput, 'Image from URL');
      onAnalysisComplete(result);
      setLinkInput('');
    } catch (err) {
      console.error('Link analysis failed:', err);
      setError(err.message || 'Failed to analyze image from URL');
    } finally {
      setIsProcessing(false);
    }
  };

  const styles = {
    container: {
      marginTop: '20px'
    },
    form: {
      display: 'flex',
      gap: '10px'
    },
    input: {
      flex: 1,
      background: 'rgba(255, 255, 255, 0.05)',
      border: '1px solid rgba(255, 255, 255, 0.2)',
      borderRadius: '4px',
      padding: '12px',
      color: '#ffffff',
      fontSize: '14px',
      outline: 'none'
    },
    button: {
      background: 'rgba(255, 255, 255, 0.1)',
      border: '1px solid rgba(255, 255, 255, 0.2)',
      borderRadius: '4px',
      padding: '12px 24px',
      color: '#ffffff',
      cursor: 'pointer',
      fontSize: '14px',
      transition: 'all 0.2s'
    },
    error: {
      color: '#ff6b6b',
      fontSize: '13px',
      marginTop: '8px'
    }
  };

  return (
    <div style={styles.container}>
      <form onSubmit={handleLinkSubmit} style={styles.form}>
        <input
          type="url"
          value={linkInput}
          onChange={(e) => setLinkInput(e.target.value)}
          placeholder="Enter image URL..."
          style={styles.input}
        />
        <button 
          type="submit" 
          style={styles.button}
        >
          Analyze
        </button>
      </form>
      {error && <div style={styles.error}>{error}</div>}
    </div>
  );
};

export default LinkSection;
