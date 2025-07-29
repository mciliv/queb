class ErrorHandler {
  static handleAIError(error, context = '') {
    console.error(`AI Error in ${context}:`, error);
    
    let errorMessage = `AI analysis failed: ${error.message}`;
    let statusCode = 500;
    
    if (error.code === 'ENOTFOUND' || error.code === 'ECONNREFUSED') {
      errorMessage = 'Network connection failed. Please check your internet connection and try again.';
      statusCode = 503;
    } else if (error.name === 'APIConnectionError') {
      errorMessage = 'Unable to connect to AI service. Please check your internet connection.';
      statusCode = 503;
    } else if (error.status === 401) {
      errorMessage = 'API authentication failed. Please check your API key configuration.';
      statusCode = 401;
    } else if (error.status === 429) {
      errorMessage = 'Rate limit exceeded. Please wait a moment and try again.';
      statusCode = 429;
    } else if (error.status === 503) {
      errorMessage = 'AI service temporarily unavailable. Please try again in a few moments.';
      statusCode = 503;
    } else if (error.message.includes('timeout')) {
      errorMessage = 'Request timeout: The AI service is taking too long to respond.';
      statusCode = 408;
    }
    
    return { errorMessage, statusCode };
  }
  
  static handleValidationError(error, context = '') {
    console.error(`Validation Error in ${context}:`, error);
    return {
      errorMessage: `Invalid input data: ${error.message}`,
      statusCode: 400
    };
  }
  
  static handleProcessingError(error, context = '') {
    console.error(`Processing Error in ${context}:`, error);
    return {
      errorMessage: `Processing failed: ${error.message}`,
      statusCode: 500
    };
  }
}

module.exports = ErrorHandler; 