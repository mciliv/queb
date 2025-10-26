// Simple logging utility to reduce noise in development
//
// Usage:
//   import logger from './logger.js';
//   logger.info('message');    // Shows in development only
//   logger.warn('warning');    // Always shows
//   logger.error('error');     // Always shows
//   logger.debug('debug');     // Shows in debug mode only
//
// Browser console controls:
//   window.debugLogger.enable()   - Enable debug logging
//   window.debugLogger.disable()  - Disable debug logging
//   window.debugLogger.verbose()  - Enable verbose (show all repetitive messages)
//   window.debugLogger.quiet()    - Disable verbose mode
//   window.debugLogger.clear()    - Clear duplicate message cache
//   window.debugLogger.status()   - Show current settings
//
// URL parameters:
//   ?debug=true     - Enable debug logging
//   ?verbose=true   - Enable verbose logging
const host = window.location.hostname;
const isLocalHost = host === 'localhost' || host === '127.0.0.1' || host === '::1' || /\.local$/i.test(host);
const isDebugMode = isLocalHost || 
                   window.location.search.includes('debug=true') ||
                   localStorage.getItem('debug-logging') === 'true';

const isVerboseMode = window.location.search.includes('verbose=true') ||
                     localStorage.getItem('verbose-logging') === 'true';

// Track logged messages to avoid spam
const loggedMessages = new Map();
const DUPLICATE_THRESHOLD = 3; // Allow up to 3 duplicates before suppressing

// Send logs to backend
const sendToBackend = async (level, message, source) => {
  // Only send errors and warnings to backend
  if (level !== 'error' && level !== 'warn') {
    return;
  }

  try {
    await fetch('/api/log-error', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        type: level,
        message: String(message),
        source: source || 'client',
        timestamp: new Date().toISOString()
      })
    });
  } catch (e) {
    // Silently fail - don't create error loops
  }
};

const createLogger = (level, prefix, consoleFn) => {
  return (message, ...args) => {
    // In production, only show errors and warnings
    if (!isDebugMode && (level === 'info' || level === 'log')) {
      return;
    }

    // For verbose logging, check if we've logged this message too many times
    if (!isVerboseMode && level === 'info') {
      const messageKey = typeof message === 'string' ? message : String(message);
      const count = (loggedMessages.get?.(messageKey) || 0) + 1;

      if (count > DUPLICATE_THRESHOLD) {
        return; // Suppress repetitive messages
      }

      loggedMessages.set(messageKey, count);
    }

    if (prefix) {
      consoleFn(`${prefix} ${message}`, ...args);
    } else {
      consoleFn(message, ...args);
    }

    // Send errors and warnings to backend
    sendToBackend(level, message, args[0]?.source);
  };
};

const logger = {
  error: createLogger('error', '❌', console.error),
  warn: createLogger('warn', '⚠️', console.warn),
  info: createLogger('info', 'ℹ️', console.info),
  log: createLogger('log', '', console.log),
  debug: createLogger('debug', '🐛', console.log),
  
  // Utility methods
  setDebug: (enabled) => {
    localStorage.setItem('debug-logging', enabled ? 'true' : 'false');
  },
  
  setVerbose: (enabled) => {
    localStorage.setItem('verbose-logging', enabled ? 'true' : 'false');
  },
  
  clearMessageCache: () => {
    loggedMessages.clear();
  }
};

// Add global debug helpers (silent in production)
if (isDebugMode) {
  window.debugLogger = {
    enable: () => logger.setDebug(true),
    disable: () => logger.setDebug(false),
    verbose: () => logger.setVerbose(true),
    quiet: () => logger.setVerbose(false),
    clear: () => logger.clearMessageCache(),
    status: () => console.log({
      debug: localStorage.getItem('debug-logging'),
      verbose: localStorage.getItem('verbose-logging'),
      isDev: isDebugMode
    })
  };
}

export default logger;
