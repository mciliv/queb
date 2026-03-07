require('dotenv').config();
const AIService = require('./AIService');

async function testAgents() {
  const ai = new AIService();
  
  console.log('--- Testing Specialized Agents ---');
  
  try {
    console.log('Running as Coder...');
    const coderResp = await ai.runAgentTask('coder', 'Write a hello world function in Python.');
    console.log('Coder response:', coderResp.content || coderResp);
    
    console.log('
Running as Architect...');
    const archResp = await ai.runAgentTask('architect', 'Suggest a cloud architecture for a global real-time chat app.');
    console.log('Architect response:', archResp.content || archResp);
    
    console.log('
--- Testing Failover (Simulation) ---');
    // We can simulate a 429 by mocking callAPI once
    const originalCallAPI = ai.callAPI;
    let failCount = 0;
    
    ai.callAPI = async (params) => {
      if (failCount < 1) {
        failCount++;
        console.log('Simulating 429 Rate Limit for current provider...');
        throw new Error('429 Too Many Requests');
      }
      return originalCallAPI.call(ai, params);
    };
    
    console.log('Testing failover from OpenAI to xAI...');
    const failoverResp = await ai.runAgentTask('coder', 'Tell me a joke about robots.');
    console.log('Failover successful. Response from backup:', failoverResp.content || failoverResp);
    
  } catch (error) {
    console.error('Test failed:', error.message);
  }
}

testAgents();
