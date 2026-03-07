#!/usr/bin/env node
/**
 * AI CLI Tool
 * Usage: node i/ask.js [--agent <coder|planner|fast>] "Your question here"
 */

require('dotenv').config();
const { aiManager } = require('./AIManager.cjs');
const { report } = require('./lib/apm/client');

async function main() {
  const args = process.argv.slice(2);

  if (args.length === 0 || args.includes('--help') || args.includes('-h')) {
    console.log('Usage: node i/ask.js [--agent <type>] "Your prompt"');
    console.log('\nAgents:');
    console.log('  coder    - High-fidelity code (Priority: OpenAI -> xAI -> Google)');
    console.log('  planner  - Deep reasoning/Context (Priority: Google -> OpenAI -> xAI)');
    console.log('  fast     - Quick responses (Priority: Google -> OpenAI)');
    process.exit(0);
  }

  let agent = 'fast';
  let promptIndex = 0;

  if (args[0] === '--agent' || args[0] === '-a') {
    agent = args[1];
    promptIndex = 2;
  }

  const prompt = args.slice(promptIndex).join(' ');

  if (!prompt) {
    console.error('Error: No prompt provided.');
    process.exit(1);
  }

  const start = Date.now();

  try {
    report({ name: 'ask', status: 'running', agent, prompt: prompt.slice(0, 200) });
    const result = await aiManager.execute(agent, { prompt });

    console.log('\n' + '-'.repeat(40));
    console.log(result.text);
    console.log('-'.repeat(40));
    console.log(`Responded via: ${result.providerUsed} (${result.modelUsed})\n`);

    report({ name: 'ask', status: 'idle', agent, provider: result.providerUsed, model: result.modelUsed, durationMs: Date.now() - start });
  } catch (error) {
    report({ name: 'ask', status: 'error', agent, error: error.message });
    console.error(`\nError: ${error.message}`);
    process.exit(1);
  }
}

main();
