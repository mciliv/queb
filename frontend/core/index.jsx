// Forward browser console logs to backend so they appear in Cursor terminal
if (process.env.NODE_ENV !== 'production') {
  const forward = (type, args) => {
    try {
      fetch('/api/log-error', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          type,
          message: args.map(a => (typeof a === 'string' ? a : JSON.stringify(a))).join(' '),
          timestamp: new Date().toISOString(),
          source: 'frontend'
        })
      });
    } catch (e) {}
  };

  const originalLog = console.log;
  const originalWarn = console.warn;
  const originalError = console.error;

  console.log = (...args) => {
    forward('log', args);
    originalLog(...args);
  };
  console.warn = (...args) => {
    forward('warn', args);
    originalWarn(...args);
  };
  console.error = (...args) => {
    forward('error', args);
    originalError(...args);
  };
}

import React from 'react';
import ReactDOM from 'react-dom/client';
import App from './App.jsx';

const root = ReactDOM.createRoot(document.getElementById('root'));
root.render(
  <React.StrictMode>
    <App />
  </React.StrictMode>
);