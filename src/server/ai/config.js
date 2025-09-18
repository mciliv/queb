// Central AI configuration — single provider, single model
const path = require('path');
require('dotenv').config({
  path: path.resolve(process.cwd(), '.env'),
  debug: false,
  quiet: true
});

module.exports = {
  provider: 'openai',
  baseURL: process.env.OPENAI_BASE_URL || 'https://api.openai.com/v1',
  apiKey: process.env.OPENAI_API_KEY || '',
  model: process.env.OPENAI_MODEL || process.env.OPENAI_DEFAULT_MODEL || 'gpt-4o',
};


