// Central AI configuration — single provider, single model
const path = require('path');
require('dotenv').config({
  path: path.resolve(process.cwd(), '.env'),
  debug: false,
  quiet: true
});

const config = {
  provider: 'openai',
  baseURL: process.env.OPENAI_BASE_URL || 'https://api.openai.com/v1',
  apiKey: process.env.OPENAI_API_KEY || '',
  model: process.env.OPENAI_MODEL || process.env.OPENAI_DEFAULT_MODEL || 'gpt-4o',
};

// In test mode, ensure a harmless default even if env vars are missing
if (process.env.NODE_ENV === 'test') {
  if (!config.apiKey) config.apiKey = 'test-key';
  if (!config.model) config.model = 'gpt-4o';
}

module.exports = config;


