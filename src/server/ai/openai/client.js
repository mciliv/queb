// Modern OpenAI client setup (CommonJS compatible) with broad export support
let OpenAIModule = null;
try {
  OpenAIModule = require('openai');
} catch (_) {
  OpenAIModule = null;
}

const aiConfig = require('../config');

function createClient() {
  if (!OpenAIModule) return null;

  // Support various export shapes:
  // - ESM default export
  // - { OpenAI } named export
  // - legacy constructor export
  const CandidateCtor = OpenAIModule?.OpenAI || OpenAIModule?.default || OpenAIModule;

  // In tests, allow a stub object to be used directly
  if (typeof CandidateCtor !== 'function') {
    return CandidateCtor; // assume already a client-like object (mock)
  }

  const options = {
    apiKey: aiConfig.apiKey,
    timeout: 30000,
    maxRetries: 2,
  };

  if (aiConfig.baseURL) {
    options.baseURL = aiConfig.baseURL;
  }

  return new CandidateCtor(options);
}

module.exports = {
  createClient
};


