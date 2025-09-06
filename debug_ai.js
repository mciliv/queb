const aiConfig = require('./backend/ai/config');
const { createClient } = require('./backend/ai/openai/client');

console.log('AI Config:');
console.log('- API Key exists:', !!aiConfig.apiKey);
console.log('- API Key length:', aiConfig.apiKey ? aiConfig.apiKey.length : 0);
console.log('- Model:', aiConfig.model);

console.log('\nOpenAI Client:');
const client = createClient();
console.log('- Client created:', !!client);

if (client) {
  console.log('- Client type:', typeof client);
  console.log('- Client has chat property:', !!client.chat);
}

console.log('\nEnvironment check:');
console.log('- OPENAI_API_KEY env:', !!process.env.OPENAI_API_KEY);
console.log('- OPENAI_API_KEY length:', process.env.OPENAI_API_KEY ? process.env.OPENAI_API_KEY.length : 0);
