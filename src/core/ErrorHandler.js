/**
 * Unified Error Handling System
 * 
 * This module provides a deep interface for error handling, logging, and recovery.
 * It hides the complexity of error classification, logging strategies, and user-friendly messaging.
 * 
 * Philosophy: "Exception handling is one of the worst sources of complexity in software systems"
 * This module reduces that complexity by providing consistent error handling patterns.
 */

class ErrorHandler {
  constructor() {
    this.logger = null;
    this.errorCounts = new Map();
    this.lastErrors = new Map();
    this.initialized = false;
  }

  /**
   * Initialize error handler with logger
   */
  initialize(logger = console) {
    // Ensure logger has all required methods
    const bind = (obj, method, fallback) => {
      const fn = obj && typeof obj[method] === 'function' ? obj[method] : fallback;
      return typeof fn === 'function' ? fn.bind(obj) : fallback.bind(console);
    };

    this.logger = {
      log: bind(logger, 'log', console.log),
      info: bind(logger, 'info', console.info || console.log),
      warn: bind(logger, 'warn', console.warn || console.log),
      error: bind(logger, 'error', console.error || console.log),
      debug: bind(logger, 'debug', console.debug || console.log)
    };
    this.initialized = true;
    
    // Set up global error handlers in non-test environments
    const isTest = process.env.NODE_ENV === 'test' || process.env.JEST_WORKER_ID;
    if (!isTest) {
      this._setupGlobalHandlers();
    }
  }

  /**
   * Handle application errors with context and recovery
   * Deep module: complex error processing behind simple interface
   */
  async handle(error, context = {}) {
    if (!this.initialized) {
      this.initialize();
    }

    const errorInfo = this._analyzeError(error, context);
    const userMessage = this._generateUserMessage(errorInfo);
    
    // Log the error
    this._logError(errorInfo);
    
    // Track error frequency
    this._trackError(errorInfo);
    
    // Attempt recovery if possible
    const recovery = this._attemptRecovery(errorInfo);
    
    return {
      message: userMessage,
      code: errorInfo.code,
      recoverable: recovery.possible,
      recovery: recovery.action,
      timestamp: errorInfo.timestamp,
      context: context.userFacing ? errorInfo.context : undefined
    };
  }

  /**
   * Handle AI service errors specifically
   */
  handleAIError(error, context = {}) {
    const enhancedContext = {
      ...context,
      category: 'ai_service',
      service: 'openai'
    };

    const result = this.handle(error, enhancedContext);
    
    // AI-specific recovery suggestions
    if (error.message?.includes('API key')) {
      result.recovery = 'Check your OpenAI API key configuration';
    } else if (error.message?.includes('rate limit')) {
      result.recovery = 'Wait a moment and try again';
    } else if (error.message?.includes('timeout')) {
      result.recovery = 'The AI service is busy. Please try again';
    }
    
    return result;
  }

  /**
   * Handle network errors
   */
  handleNetworkError(error, context = {}) {
    const enhancedContext = {
      ...context,
      category: 'network',
      retryable: true
    };

    const result = this.handle(error, enhancedContext);
    result.recovery = 'Check your internet connection and try again';
    
    return result;
  }

  /**
   * Handle validation errors
   */
  handleValidationError(error, context = {}) {
    const enhancedContext = {
      ...context,
      category: 'validation',
      userCaused: true
    };

    return this.handle(error, enhancedContext);
  }

  /**
   * Handle file system errors
   */
  handleFileSystemError(error, context = {}) {
    const enhancedContext = {
      ...context,
      category: 'filesystem'
    };

    const result = this.handle(error, enhancedContext);
    
    if (error.code === 'ENOENT') {
      result.recovery = 'Required file not found. Please check the installation';
    } else if (error.code === 'EACCES') {
      result.recovery = 'Permission denied. Check file permissions';
    } else if (error.code === 'ENOSPC') {
      result.recovery = 'Disk full. Please free up space';
    }
    
    return result;
  }

  /**
   * PRIVATE IMPLEMENTATION - Complex error processing logic
   */

  /**
   * Analyze error and extract useful information
   */
  _analyzeError(error, context) {
    const timestamp = new Date().toISOString();
    const errorId = this._generateErrorId(error, context);
    
    // Extract error details
    let message = 'Unknown error occurred';
    let stack = null;
    let code = 'UNKNOWN_ERROR';
    
    if (error instanceof Error) {
      message = error.message;
      stack = error.stack;
      code = error.code || error.name || 'ERROR';
    } else if (typeof error === 'string') {
      message = error;
      code = 'STRING_ERROR';
    } else if (error && typeof error === 'object') {
      message = error.message || error.error || JSON.stringify(error);
      code = error.code || 'OBJECT_ERROR';
    }

    // Classify error severity
    const severity = this._classifyErrorSeverity(message, context);
    
    // Determine if error is user-facing
    const isUserFacing = this._isUserFacingError(message, context);
    
    return {
      id: errorId,
      message,
      originalError: error,
      stack,
      code,
      severity,
      timestamp,
      context: {
        ...context,
        userAgent: typeof window !== 'undefined' ? window.navigator?.userAgent : undefined,
        url: typeof window !== 'undefined' ? window.location?.href : undefined
      },
      isUserFacing,
      category: context.category || this._categorizeError(message)
    };
  }

  /**
   * Generate user-friendly error messages
   */
  _generateUserMessage(errorInfo) {
    // Use user-friendly messages for common errors
    const friendlyMessages = {
      'AI service unavailable': 'AI analysis is temporarily unavailable. Please try again later.',
      'Network error': 'Connection failed. Please check your internet connection.',
      'Invalid input': 'Please check your input and try again.',
      'File not found': 'Required file could not be found.',
      'Permission denied': 'Access denied. Please check permissions.',
      'Timeout': 'Request timed out. Please try again.',
      'Rate limit': 'Too many requests. Please wait a moment and try again.'
    };

    // Check for known error patterns
    for (const [pattern, friendlyMessage] of Object.entries(friendlyMessages)) {
      if (errorInfo.message.toLowerCase().includes(pattern.toLowerCase())) {
        return friendlyMessage;
      }
    }

    // Fallback to processed error message
    if (errorInfo.isUserFacing) {
      return this._sanitizeErrorMessage(errorInfo.message);
    }

    // Generic message for technical errors
    return 'An unexpected error occurred. Please try again.';
  }

  /**
   * Classify error severity
   */
  _classifyErrorSeverity(message, context) {
    const lowerMessage = message.toLowerCase();
    
    // Critical errors
    if (lowerMessage.includes('crash') || lowerMessage.includes('fatal') || 
        lowerMessage.includes('out of memory')) {
      return 'critical';
    }
    
    // High severity
    if (lowerMessage.includes('database') || lowerMessage.includes('connection') ||
        lowerMessage.includes('authentication')) {
      return 'high';
    }
    
    // Medium severity
    if (lowerMessage.includes('timeout') || lowerMessage.includes('rate limit') ||
        lowerMessage.includes('validation')) {
      return 'medium';
    }
    
    // Low severity for user input errors
    if (context.userCaused || lowerMessage.includes('invalid input')) {
      return 'low';
    }
    
    return 'medium';
  }

  /**
   * Categorize errors for better handling
   */
  _categorizeError(message) {
    const lowerMessage = message.toLowerCase();
    
    if (lowerMessage.includes('network') || lowerMessage.includes('fetch') ||
        lowerMessage.includes('connection')) {
      return 'network';
    }
    
    if (lowerMessage.includes('api key') || lowerMessage.includes('authentication') ||
        lowerMessage.includes('unauthorized')) {
      return 'authentication';
    }
    
    if (lowerMessage.includes('validation') || lowerMessage.includes('invalid')) {
      return 'validation';
    }
    
    if (lowerMessage.includes('file') || lowerMessage.includes('path') ||
        lowerMessage.includes('directory')) {
      return 'filesystem';
    }
    
    if (lowerMessage.includes('timeout') || lowerMessage.includes('rate limit')) {
      return 'service_limit';
    }
    
    return 'application';
  }

  /**
   * Determine if error should be shown to users
   */
  _isUserFacingError(message, context) {
    if (context.userFacing !== undefined) {
      return context.userFacing;
    }
    
    // Technical errors should not be shown directly
    const technicalPatterns = [
      'stack trace', 'undefined is not a function', 'cannot read property',
      'module not found', 'syntax error', 'reference error'
    ];
    
    const lowerMessage = message.toLowerCase();
    return !technicalPatterns.some(pattern => lowerMessage.includes(pattern));
  }

  /**
   * Attempt error recovery
   */
  _attemptRecovery(errorInfo) {
    const category = errorInfo.category;
    
    switch (category) {
      case 'network':
        return { possible: true, action: 'retry_with_backoff' };
      
      case 'service_limit':
        return { possible: true, action: 'retry_after_delay' };
      
      case 'validation':
        return { possible: true, action: 'user_correction' };
      
      case 'authentication':
        return { possible: false, action: 'check_configuration' };
      
      default:
        return { possible: false, action: 'manual_intervention' };
    }
  }

  /**
   * Log error with appropriate level
   */
  _logError(errorInfo) {
    const logData = {
      id: errorInfo.id,
      message: errorInfo.message,
      code: errorInfo.code,
      severity: errorInfo.severity,
      category: errorInfo.category,
      timestamp: errorInfo.timestamp,
      context: errorInfo.context
    };

    switch (errorInfo.severity) {
      case 'critical':
        this.logger.error('🚨 CRITICAL ERROR:', logData);
        break;
      case 'high':
        this.logger.error('❌ HIGH SEVERITY ERROR:', logData);
        break;
      case 'medium':
        this.logger.warn('⚠️ ERROR:', logData);
        break;
      case 'low':
        this.logger.info('ℹ️ Low severity error:', logData);
        break;
      default:
        this.logger.error('❌ ERROR:', logData);
    }

    // Include stack trace for development
    const isDev = process.env.NODE_ENV === 'development';
    if (isDev && errorInfo.stack) {
      this.logger.debug('Stack trace:', errorInfo.stack);
    }
  }

  /**
   * Track error frequency for monitoring
   */
  _trackError(errorInfo) {
    const key = `${errorInfo.code}:${errorInfo.category}`;
    const count = this.errorCounts.get(key) || 0;
    this.errorCounts.set(key, count + 1);
    this.lastErrors.set(key, errorInfo.timestamp);
    
    // Alert on high frequency errors
    if (count > 10) {
      this.logger.warn(`⚠️ High frequency error detected: ${key} (${count + 1} occurrences)`);
    }
  }

  /**
   * Generate unique error ID
   */
  _generateErrorId(error, context) {
    const timestamp = Date.now();
    const contextKey = context.category || 'unknown';
    const errorMessage = error?.message || error || 'unknown';
    const errorKey = typeof errorMessage === 'string' 
      ? errorMessage.substring(0, 20) 
      : String(errorMessage).substring(0, 20);
    return `${contextKey}_${errorKey.replace(/\W/g, '_')}_${timestamp}`;
  }

  /**
   * Sanitize error messages for user display
   */
  _sanitizeErrorMessage(message) {
    // Remove technical details
    return message
      .replace(/at \w+\.\w+.*$/gm, '') // Remove stack trace references
      .replace(/\b[A-Z][a-z]*Error:\s*/g, '') // Remove error type prefixes
      .replace(/\s+/g, ' ') // Normalize whitespace
      .trim();
  }

  /**
   * Setup global error handlers
   */
  _setupGlobalHandlers() {
    if (typeof window !== 'undefined') {
      // Browser environment
      window.addEventListener('error', (event) => {
        this.handle(event.error, { category: 'javascript', global: true });
      });

      window.addEventListener('unhandledrejection', (event) => {
        this.handle(event.reason, { category: 'promise', global: true });
      });
    } else if (typeof process !== 'undefined') {
      // Node.js environment
      process.on('uncaughtException', (error) => {
        const result = this.handle(error, { category: 'uncaught_exception', global: true });
        console.error('Uncaught Exception:', result);
        process.exit(1);
      });

      process.on('unhandledRejection', (reason) => {
        const result = this.handle(reason, { category: 'unhandled_rejection', global: true });
        console.error('Unhandled Rejection:', result);
      });
    }
  }

  /**
   * Get error statistics for monitoring
   */
  getErrorStats() {
    return {
      totalErrors: Array.from(this.errorCounts.values()).reduce((sum, count) => sum + count, 0),
      errorsByCategory: Object.fromEntries(this.errorCounts),
      recentErrors: Object.fromEntries(this.lastErrors)
    };
  }
}

module.exports = ErrorHandler;


