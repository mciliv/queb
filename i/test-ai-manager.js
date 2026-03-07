/**
 * Test script for AIManager
 * Demonstrates model priority, specialized agents, and failover.
 */

import { aiManager } from './AIManager.js';
import 'dotenv/config';

async function runTests() {
  console.log('--- Testing AIManager Specialized Agents ---
');

  try {
    // 1. Test Coder Agent (Priority: OpenAI -> xAI -> Google)
    console.log('Task: Write a simple React component (Agent: coder)');
    const coderResp = await aiManager.execute('coder', {
      prompt: 'Write a React component for a Todo list.'
    });
    console.log(`✅ Success via ${coderResp.providerUsed} (${coderResp.modelUsed})`);
    console.log(`Content length: ${coderResp.text.length}
`);

    // 2. Test Planner Agent (Priority: Google -> OpenAI -> xAI)
    console.log('Task: Long-term technical roadmap (Agent: planner)');
    const plannerResp = await aiManager.execute('planner', {
      prompt: 'Outline a 12-month roadmap for migrating a monolith to microservices.'
    });
    console.log(`✅ Success via ${plannerResp.providerUsed} (${plannerResp.modelUsed})`);
    console.log(`Content length: ${plannerResp.text.length}
`);

    // 3. Failover Demo (Simulation)
    console.log('--- Testing Failover Logic ---
');
    console.log('Simulating 429 for OpenAI to see it fail over to xAI/Google...');
    
    // We'll mark OpenAI as throttled manually for this test
    aiManager._markThrottled('openai');
    
    const failoverResp = await aiManager.execute('coder', {
      prompt: 'What is the fastest way to learn TypeScript?'
    });
    console.log(`✅ Failover Success via ${failoverResp.providerUsed} (${failoverResp.modelUsed})`);
    
  } catch (error) {
    console.error('❌ Test failed:', error.message);
  }
}

runTests();
