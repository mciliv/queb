import React from 'react';

// Enhanced loading spinner with different states
export const LoadingSpinner = ({ size = 'medium', message = '', progress = null }) => {
  const sizeClasses = {
    small: 'w-4 h-4',
    medium: 'w-6 h-6', 
    large: 'w-8 h-8'
  };

  return (
    <div className="loading-spinner-container" style={{ textAlign: 'center', padding: '1rem' }}>
      <div 
        className={`loading-spinner ${sizeClasses[size]}`}
        style={{
          width: size === 'small' ? '16px' : size === 'large' ? '32px' : '24px',
          height: size === 'small' ? '16px' : size === 'large' ? '32px' : '24px',
          border: '2px solid #f3f3f3',
          borderTop: '2px solid #3498db',
          borderRadius: '50%',
          animation: 'spin 1s linear infinite',
          margin: '0 auto'
        }}
        role="status"
        aria-label="Loading"
      />
      
      {progress !== null && (
        <div className="progress-bar" style={{ marginTop: '0.5rem' }}>
          <div 
            className="progress-fill"
            style={{
              width: '100%',
              height: '4px',
              backgroundColor: '#f3f3f3',
              borderRadius: '2px',
              overflow: 'hidden'
            }}
          >
            <div
              style={{
                width: `${Math.min(100, Math.max(0, progress))}%`,
                height: '100%',
                backgroundColor: '#3498db',
                transition: 'width 0.3s ease',
                borderRadius: '2px'
              }}
            />
          </div>
          <div style={{ fontSize: '0.75rem', marginTop: '0.25rem', color: '#666' }}>
            {Math.round(progress)}%
          </div>
        </div>
      )}
      
      {message && (
        <div 
          className="loading-message" 
          style={{ 
            marginTop: '0.5rem', 
            fontSize: '0.875rem', 
            color: '#666',
            fontStyle: 'italic'
          }}
        >
          {message}
        </div>
      )}
    </div>
  );
};

// Skeleton loader for cards/results
export const SkeletonLoader = ({ lines = 3, height = 'auto' }) => {
  return (
    <div 
      className="skeleton-loader"
      style={{ 
        padding: '1rem',
        height,
        display: 'flex',
        flexDirection: 'column',
        gap: '0.5rem'
      }}
    >
      {Array(lines).fill(0).map((_, i) => (
        <div
          key={i}
          style={{
            height: '1rem',
            backgroundColor: '#f3f3f3',
            borderRadius: '4px',
            animation: 'pulse 1.5s ease-in-out infinite',
            width: i === lines - 1 ? '60%' : '100%'
          }}
        />
      ))}
    </div>
  );
};

// Processing status indicator
export const ProcessingStatus = ({ stage, stages = [], currentProgress = 0 }) => {
  const defaultStages = [
    'Analyzing input',
    'Converting to SMILES',
    'Generating 3D structures',
    'Rendering molecules'
  ];
  
  const stageList = stages.length > 0 ? stages : defaultStages;
  const currentStageIndex = stageList.indexOf(stage);

  return (
    <div className="processing-status" style={{ padding: '1rem' }}>
      <div style={{ marginBottom: '1rem' }}>
        <strong>Processing...</strong>
      </div>
      
      <div className="stage-list">
        {stageList.map((stageName, index) => {
          const isActive = index === currentStageIndex;
          const isCompleted = index < currentStageIndex;
          
          return (
            <div
              key={stageName}
              className={`stage-item ${isActive ? 'active' : ''} ${isCompleted ? 'completed' : ''}`}
              style={{
                display: 'flex',
                alignItems: 'center',
                padding: '0.5rem 0',
                color: isCompleted ? '#28a745' : isActive ? '#007bff' : '#6c757d'
              }}
            >
              <div
                className="stage-icon"
                style={{
                  width: '20px',
                  height: '20px',
                  borderRadius: '50%',
                  marginRight: '0.75rem',
                  display: 'flex',
                  alignItems: 'center',
                  justifyContent: 'center',
                  fontSize: '0.75rem',
                  backgroundColor: isCompleted ? '#28a745' : isActive ? '#007bff' : '#e9ecef',
                  color: isCompleted || isActive ? 'white' : '#6c757d'
                }}
              >
                {isCompleted ? '✓' : isActive ? '⋯' : index + 1}
              </div>
              <span style={{ fontSize: '0.875rem' }}>{stageName}</span>
              {isActive && (
                <LoadingSpinner size="small" />
              )}
            </div>
          );
        })}
      </div>
      
      {currentProgress > 0 && (
        <div style={{ marginTop: '1rem' }}>
          <div 
            style={{
              width: '100%',
              height: '6px',
              backgroundColor: '#e9ecef',
              borderRadius: '3px',
              overflow: 'hidden'
            }}
          >
            <div
              style={{
                width: `${currentProgress}%`,
                height: '100%',
                backgroundColor: '#007bff',
                transition: 'width 0.3s ease',
                borderRadius: '3px'
              }}
            />
          </div>
        </div>
      )}
    </div>
  );
};

// Error state with retry functionality
export const ErrorState = ({ 
  error, 
  onRetry, 
  retryCount = 0, 
  maxRetries = 3,
  showDetails = false 
}) => {
  const canRetry = retryCount < maxRetries && typeof onRetry === 'function';

  return (
    <div 
      className="error-state"
      style={{
        padding: '1rem',
        border: '1px solid #dc3545',
        borderRadius: '4px',
        backgroundColor: '#f8d7da',
        color: '#721c24',
        textAlign: 'center'
      }}
    >
      <div style={{ fontSize: '1.25rem', marginBottom: '0.5rem' }}>⚠️</div>
      <div style={{ fontWeight: 'bold', marginBottom: '0.5rem' }}>
        Something went wrong
      </div>
      <div style={{ fontSize: '0.875rem', marginBottom: '1rem' }}>
        {error?.message || error || 'An unexpected error occurred'}
      </div>
      
      {canRetry && (
        <button
          onClick={onRetry}
          style={{
            backgroundColor: '#dc3545',
            color: 'white',
            border: 'none',
            padding: '0.5rem 1rem',
            borderRadius: '4px',
            cursor: 'pointer',
            fontSize: '0.875rem',
            marginRight: '0.5rem'
          }}
        >
          Retry ({maxRetries - retryCount} attempts left)
        </button>
      )}
      
      {showDetails && error?.stack && (
        <details style={{ marginTop: '1rem', textAlign: 'left' }}>
          <summary style={{ cursor: 'pointer', fontSize: '0.75rem' }}>
            Show technical details
          </summary>
          <pre style={{ 
            fontSize: '0.625rem', 
            backgroundColor: '#fff', 
            padding: '0.5rem',
            borderRadius: '2px',
            overflow: 'auto',
            marginTop: '0.5rem'
          }}>
            {error.stack}
          </pre>
        </details>
      )}
    </div>
  );
};

// Success state
export const SuccessState = ({ message = 'Completed successfully!', children }) => {
  return (
    <div
      className="success-state"
      style={{
        padding: '1rem',
        border: '1px solid #28a745',
        borderRadius: '4px',
        backgroundColor: '#d4edda',
        color: '#155724',
        textAlign: 'center'
      }}
    >
      <div style={{ fontSize: '1.25rem', marginBottom: '0.5rem' }}>✅</div>
      <div style={{ fontWeight: 'bold', marginBottom: '0.5rem' }}>
        {message}
      </div>
      {children}
    </div>
  );
};

// Add CSS animations via style injection
if (typeof document !== 'undefined') {
  const style = document.createElement('style');
  style.textContent = `
    @keyframes spin {
      0% { transform: rotate(0deg); }
      100% { transform: rotate(360deg); }
    }
    
    @keyframes pulse {
      0%, 100% { opacity: 1; }
      50% { opacity: 0.5; }
    }
    
    .loading-spinner-container .loading-spinner {
      animation: spin 1s linear infinite;
    }
    
    .skeleton-loader > div {
      animation: pulse 1.5s ease-in-out infinite;
    }
  `;
  document.head.appendChild(style);
}