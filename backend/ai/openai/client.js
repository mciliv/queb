let OpenAIClient = null;
try { OpenAIClient = require('openai').OpenAI; } catch (_) { OpenAIClient = null; }
const aiConfig = require('../config');

function createClient() {
  if (!OpenAIClient) return null;
  const options = { apiKey: aiConfig.apiKey };
  if (aiConfig.baseURL) options.baseURL = aiConfig.baseURL;
  return new OpenAIClient(options);
}

module.exports = {
  createClient,
};


