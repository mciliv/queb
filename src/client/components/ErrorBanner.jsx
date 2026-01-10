import React from 'react';

const ErrorBanner = ({ error, onDismiss }) => {
  if (!error) return null;
  
  let errorMessage;
  if (typeof error === 'string') {
    errorMessage = error;
  } else if (error?.message) {
    errorMessage = error.message;
  } else if (error?.code) {
    errorMessage = `Error ${error.code}: ${JSON.stringify(error)}`;
  } else {
    errorMessage = `Error: ${JSON.stringify(error)}`;
  }
  
  return (
    <div className="error-banner">
      <div className="error-banner-content">
        <div className="error-banner-header">
          <span className="error-icon">⚠️</span>
          <div className="error-message">{errorMessage}</div>
          {onDismiss && (
            <button
              onClick={onDismiss}
              className="error-dismiss"
              aria-label="Dismiss error"
            >
              ×
            </button>
          )}
        </div>
      </div>
    </div>
  );
};

export default ErrorBanner;