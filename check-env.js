// Quick environment checker
require('dotenv').config();

console.log('🔍 Environment Check:');
console.log('  NODE_ENV:', process.env.NODE_ENV);
console.log('  PORT:', process.env.PORT);
console.log('  OPENAI_API_KEY present:', !!process.env.OPENAI_API_KEY);
console.log('  OPENAI_API_KEY length:', process.env.OPENAI_API_KEY ? process.env.OPENAI_API_KEY.length : 0);

if (!process.env.OPENAI_API_KEY) {
  console.log('❌ OPENAI_API_KEY not found in environment');
  console.log('💡 Check if .env file exists and contains OPENAI_API_KEY=your_key_here');
} else {
  console.log('✅ OPENAI_API_KEY loaded successfully');
}