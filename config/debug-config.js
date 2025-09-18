// Debug config loading
const config = require('./env.js');

console.log('🔧 Config Debug:');
console.log('  config.NODE_ENV:', config.NODE_ENV);
console.log('  config.OPENAI_API_KEY present:', !!config.OPENAI_API_KEY);
console.log('  config.OPENAI_API_KEY length:', config.OPENAI_API_KEY ? config.OPENAI_API_KEY.length : 0);
console.log('  process.env.OPENAI_API_KEY present:', !!process.env.OPENAI_API_KEY);

// Test the server logic
const openaiApiKey = config.OPENAI_API_KEY;
console.log('  Final openaiApiKey present:', !!openaiApiKey);