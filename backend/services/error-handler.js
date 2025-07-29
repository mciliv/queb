/**
 * CENTRALIZED ERROR HANDLER
 * Purpose: Unified error handling for all backend operations
 * Provides: Consistent error messages, status codes, and logging
 */

class ErrorHandler {
  /**
   * Handle AI service errors with consistent messaging
   * @param {Error} error - The original error
   * @param {string} context - Where the error occurred
   * @returns {Object} Standardized error response
   */
  static handleAIError(error, context = '') {
    console.error(`🚨 AI Error in ${context}:`, error);
    
    let errorMessage = `AI analysis failed: ${error.message}`;
    let statusCode = 500;
    
    // Network-related errors
    if (error.code === 'ENOTFOUND' || error.code === 'ECONNREFUSED') {
      errorMessage = 'Network connection failed. Please check your internet connection and try again.';
      statusCode = 503;
    } else if (error.name === 'APIConnectionError') {
      errorMessage = 'Unable to connect to AI service. Please check your internet connection.';
      statusCode = 503;
    } 
    // API-specific errors
    else if (error.status === 401) {
      errorMessage = 'API authentication failed. Please check your API key configuration.';
      statusCode = 401;
    } else if (error.status === 429) {
      errorMessage = 'Rate limit exceeded. Please wait a moment and try again.';
      statusCode = 429;
    } else if (error.status === 503) {
      errorMessage = 'AI service temporarily unavailable. Please try again in a few moments.';
      statusCode = 503;
    } 
    // Timeout errors
    else if (error.message.includes('timeout')) {
      errorMessage = 'Request timeout: The AI service is taking too long to respond.';
      statusCode = 408;
    }
    
    return this.createErrorResponse(errorMessage, statusCode, context);
  }
  
  /**
   * Handle validation errors with consistent messaging
   * @param {Error} error - The validation error
   * @param {string} context - Where the validation failed
   * @returns {Object} Standardized error response
   */
  static handleValidationError(error, context = '') {
    console.error(`🚨 Validation Error in ${context}:`, error);
    return this.createErrorResponse(
      `Invalid input data: ${error.message}`,
      400,
      context
    );
  }
  
  /**
   * Handle processing errors with consistent messaging
   * @param {Error} error - The processing error
   * @param {string} context - Where the processing failed
   * @returns {Object} Standardized error response
   */
  static handleProcessingError(error, context = '') {
    console.error(`🚨 Processing Error in ${context}:`, error);
    return this.createErrorResponse(
      `Processing failed: ${error.message}`,
      500,
      context
    );
  }

  /**
   * Handle service initialization errors
   * @param {Error} error - The initialization error
   * @param {string} serviceName - Name of the service that failed
   * @returns {Object} Standardized error response
   */
  static handleServiceError(error, serviceName = '') {
    console.error(`🚨 Service Error in ${serviceName}:`, error);
    return this.createErrorResponse(
      `Service initialization failed: ${serviceName} - ${error.message}`,
      503,
      serviceName
    );
  }

  /**
   * Create standardized error response
   * @param {string} errorMessage - Human-readable error message
   * @param {number} statusCode - HTTP status code
   * @param {string} context - Context where error occurred
   * @returns {Object} Standardized error response
   */
  static createErrorResponse(errorMessage, statusCode, context) {
    return {
      errorMessage,
      statusCode,
      context,
      timestamp: new Date().toISOString()
    };
  }

  /**
   * Log successful operations for debugging
   * @param {string} operation - Description of the operation
   * @param {string} context - Context of the operation
   * @param {Object} data - Optional data to log
   */
  static logSuccess(operation, context = '', data = null) {
    const logData = {
      operation,
      context,
      timestamp: new Date().toISOString(),
      ...(data && { data })
    };
    console.log(`✅ Success: ${operation} in ${context}`, logData);
  }
}

module.exports = ErrorHandler; 