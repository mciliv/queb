// Centralized logging system for production
class Logger {
  constructor() {
    this.logLevel = this.getLogLevel();
    this.enableConsole = this.logLevel === 'debug';
    this.serverEndpoint = '/api/logs';
    this.batchSize = 10;
    this.batchTimeout = 5000; // 5 seconds
    this.logQueue = [];
    this.batchTimer = null;
    
    // Log initialization
    this.info('Logger initialized', { 
      level: this.logLevel, 
      enableConsole: this.enableConsole,
      serverEndpoint: this.serverEndpoint 
    });
  }

  getLogLevel() {
    // Check URL params first, then localStorage, then default to 'info'
    const urlParams = new URLSearchParams(window.location.search);
    const logLevel = urlParams.get('log') || localStorage.getItem('logLevel') || 'info';
    return ['debug', 'info', 'warn', 'error'].includes(logLevel) ? logLevel : 'info';
  }

  // Set log level dynamically
  setLogLevel(level) {
    if (['debug', 'info', 'warn', 'error'].includes(level)) {
      this.logLevel = level;
      this.enableConsole = level === 'debug';
      localStorage.setItem('logLevel', level);
      
      // Log the change
      this.info('Log level changed', { newLevel: level });
    }
  }

  // Send logs to server
  async sendToServer(logs) {
    try {
      await fetch(this.serverEndpoint, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          logs,
          timestamp: new Date().toISOString(),
          userAgent: navigator.userAgent,
          url: window.location.href
        })
      });
    } catch (error) {
      // Fallback to console if server logging fails
      if (this.enableConsole) {
        console.error('Failed to send logs to server:', error);
      }
    }
  }

  // Process log queue
  processQueue() {
    if (this.logQueue.length === 0) return;

    const logs = [...this.logQueue];
    this.logQueue = [];
    
    this.sendToServer(logs);
  }

  // Schedule batch processing
  scheduleBatch() {
    if (this.batchTimer) clearTimeout(this.batchTimer);
    
    this.batchTimer = setTimeout(() => {
      this.processQueue();
    }, this.batchTimeout);
  }

  // Add log to queue
  addToQueue(level, message, data = null) {
    const logEntry = {
      level,
      message,
      data,
      timestamp: new Date().toISOString(),
      source: this.getCallerInfo()
    };

    this.logQueue.push(logEntry);

    // Process immediately if batch is full
    if (this.logQueue.length >= this.batchSize) {
      this.processQueue();
    } else {
      this.scheduleBatch();
    }

    // Also log to console in debug mode
    if (this.enableConsole) {
      const consoleMethod = level === 'error' ? 'error' : 
                           level === 'warn' ? 'warn' : 'log';
      console[consoleMethod](`[${level.toUpperCase()}]`, message, data || '');
    }
  }

  // Get caller information
  getCallerInfo() {
    try {
      const stack = new Error().stack;
      const lines = stack.split('\n');
      // Skip first 3 lines (Error, Logger method, calling method)
      const callerLine = lines[3] || '';
      const match = callerLine.match(/at\s+(.+?)\s+\((.+?):(\d+):(\d+)\)/);
      if (match) {
        return {
          function: match[1],
          file: match[2].split('/').pop(),
          line: match[3],
          column: match[4]
        };
      }
    } catch (e) {
      // Ignore errors in caller detection
    }
    return { function: 'unknown', file: 'unknown', line: '0', column: '0' };
  }

  // Public logging methods
  debug(message, data = null) {
    if (this.logLevel === 'debug') {
      this.addToQueue('debug', message, data);
    }
  }

  info(message, data = null) {
    if (['debug', 'info'].includes(this.logLevel)) {
      this.addToQueue('info', message, data);
    }
  }

  warn(message, data = null) {
    if (['debug', 'info', 'warn'].includes(this.logLevel)) {
      this.addToQueue('warn', message, data);
    }
  }

  error(message, data = null) {
    // Always log errors
    this.addToQueue('error', message, data);
  }

  // Log API calls
  apiCall(url, method, requestData = null) {
    this.info('API Call', { url, method, requestData });
  }

  apiResponse(url, status, responseData = null) {
    this.info('API Response', { url, status, responseData });
  }

  apiError(url, error) {
    this.error('API Error', { url, error: error.message, stack: error.stack });
  }

  // Log user interactions
  userAction(action, details = null) {
    this.info('User Action', { action, details });
  }

  // Log payment events
  paymentEvent(event, details = null) {
    this.info('Payment Event', { event, details });
  }

  // Log camera events
  cameraEvent(event, details = null) {
    this.info('Camera Event', { event, details });
  }

  // Log analysis events
  analysisEvent(event, details = null) {
    this.info('Analysis Event', { event, details });
  }

  // Flush remaining logs (call before page unload)
  flush() {
    this.processQueue();
  }
}

// Create global logger instance
const logger = new Logger();

// Override console methods in production
if (logger.logLevel !== 'debug') {
  // Store original console methods
  const originalConsole = {
    log: console.log,
    info: console.info,
    warn: console.warn,
    error: console.error,
    debug: console.debug
  };

  // Override console methods to use our logger
  console.log = (...args) => logger.info(args.join(' '));
  console.info = (...args) => logger.info(args.join(' '));
  console.warn = (...args) => logger.warn(args.join(' '));
  console.error = (...args) => logger.error(args.join(' '));
  console.debug = (...args) => logger.debug(args.join(' '));

  // Restore console for debugging if needed
  window.restoreConsole = () => {
    Object.assign(console, originalConsole);
  };
}

// Flush logs before page unload
window.addEventListener('beforeunload', () => {
  logger.flush();
});

// Add global logger controls
window.logger = logger;
window.setLogLevel = (level) => logger.setLogLevel(level);

// Export logger
export { logger };
export default logger; 