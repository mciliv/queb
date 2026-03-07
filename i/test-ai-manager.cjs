/**
 * Test script for AIManager (CommonJS)
 */

require('dotenv').config();
const { aiManager } = require('./AIManager.cjs');

async function runTests() {
  console.log('--- Testing AIManager (CJS) ---\n');

  try {
    // 1. Test Coder Agent
    console.log('Task: Write a simple Python script (Agent: coder)');
    const coderResp = await aiManager.execute('coder', {
      prompt: 'Write a python script that prints primes up to 100.'
    });
    console.log(`✅ Success via ${coderResp.providerUsed} (${coderResp.modelUsed})\n`);

    // 2. Test Planner Agent
    console.log('Task: Long-term technical roadmap (Agent: planner)');
    const plannerResp = await aiManager.execute('planner', {
      prompt: 'Suggest a strategy for scaling a Node.js app to 1M users.'
    });
    console.log(`✅ Success via ${plannerResp.providerUsed} (${plannerResp.modelUsed})\n`);

    // 3. Failover Demo (Simulation)
    console.log('--- Testing Failover Logic ---\n');
    console.log('Simulating 429 for "gpt-4o" to see it fail over to the next best model...');
    
    aiManager._markThrottled('gpt-4o');
    
    const failoverResp = await aiManager.execute('coder', {
      prompt: 'Why is the sky blue?'
    });
    console.log(`✅ Failover Success via ${failoverResp.providerUsed} (${failoverResp.modelUsed})`);
    
  } catch (error) {
    console.error('❌ Test failed:', error.message);
  }
}

runTests();
