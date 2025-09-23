// Modern OpenAI client setup (CommonJS compatible)
let OpenAI = null;
try { 
  OpenAI = require('openai').default || require('openai'); 
} catch (_) { 
  OpenAI = null; 
}

const aiConfig = require('../config');

function createClient() {
  if (!OpenAI) return null;
  
  const options = { 
    apiKey: aiConfig.apiKey,
    timeout: 30000, // 30 second timeout
    maxRetries: 2   // Retry failed requests twice
  };
  
  if (aiConfig.baseURL) {
    options.baseURL = aiConfig.baseURL;
  }
  
  return new OpenAI(options);
}

module.exports = {
  createClient,
  OpenAI
};


