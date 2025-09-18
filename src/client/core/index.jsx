// Forward browser console logs (with stack/location) to backend in dev
if (process.env.NODE_ENV !== 'production') {
  const parseStackTop = (stack) => {
    if (!stack || typeof stack !== 'string') return null;
    const lines = stack.split('\n').map(s => s.trim());
    // Skip the first line ("Error") and our own wrapper frames
    for (let i = 1; i < lines.length; i++) {
      const ln = lines[i];
      const m = ln.match(/\((.*):(\d+):(\d+)\)$/) || ln.match(/at (.*):(\d+):(\d+)/);
      if (m) {
        return { file: m[1], line: Number(m[2]), column: Number(m[3]) };
      }
    }
    return null;
  };

  const toMessage = (args) => args.map(a => {
    if (a instanceof Error) return a.message;
    if (typeof a === 'string') return a;
    try { return JSON.stringify(a); } catch (_) { return String(a); }
  }).join(' ');

  const capture = (type, args) => {
    try {
      let stack = '';
      // Prefer explicit Error arg stack
      const errArg = args.find(a => a instanceof Error);
      if (errArg && errArg.stack) {
        stack = errArg.stack;
      } else {
        try { throw new Error(); } catch (e) { stack = e.stack || ''; }
      }
      const loc = parseStackTop(stack);
      fetch('/api/log-error', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          type,
          message: toMessage(args),
          timestamp: new Date().toISOString(),
          source: 'frontend',
          stack,
          location: loc,
          url: (typeof window !== 'undefined' && window.location && window.location.href) ? window.location.href : ''
        })
      });
    } catch (_) {}
  };

  const originalLog = console.log;
  const originalWarn = console.warn;
  const originalError = console.error;

  console.log = (...args) => { capture('log', args); originalLog(...args); };
  console.warn = (...args) => { capture('warn', args); originalWarn(...args); };
  console.error = (...args) => { capture('error', args); originalError(...args); };

  // Also forward global errors/rejections (but suppress Chrome extension errors)
  window.addEventListener('error', (e) => {
    // Suppress Chrome extension runtime errors
    if (e.message && (e.message.includes('Extension context invalidated') ||
        e.message.includes('message port closed') ||
        e.message.includes('runtime.lastError'))) {
      return; // Don't forward Chrome extension errors
    }
    capture('error', [e.message, e.error || '']);
  });
  window.addEventListener('unhandledrejection', (e) => {
    // Suppress Chrome extension promise rejections
    if (e.reason && typeof e.reason === 'string' &&
        (e.reason.includes('Extension context invalidated') ||
         e.reason.includes('message port closed') ||
         e.reason.includes('runtime.lastError'))) {
      return; // Don't forward Chrome extension errors
    }
    capture('error', ['UnhandledRejection', e.reason || '']);
  });
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