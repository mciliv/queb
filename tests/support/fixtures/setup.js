// test/fixtures/setup.js - Jest setup file for mocking dependencies

const cleanup = require('./cleanup-registry');

// Mock OpenAI before any modules are imported with chat.completions API
jest.mock('openai', () => {
  return {
    default: jest.fn().mockImplementation(() => ({
      chat: {
        completions: {
          create: jest.fn().mockImplementation(async ({ messages }) => {
            const last = Array.isArray(messages) ? [...messages].reverse().find(m => m.role === 'user') : null;
            const content = last?.content;
            let text = '';
            if (typeof content === 'string') {
              text = content;
            } else if (Array.isArray(content)) {
              const textPart = content.find(p => p && p.type === 'text');
              text = textPart?.text || '';
            }

            const lower = (text || '').toLowerCase();
            const chemicals = [];
            if (lower.includes('water')) chemicals.push({ name: 'water', smiles: 'O' });
            if (lower.includes('ethanol') || lower.includes('wine')) chemicals.push({ name: 'ethanol', smiles: 'CCO' });
            if (lower.includes('sodium') && lower.includes('chloride')) chemicals.push({ name: 'sodium chloride', smiles: '[Na+].[Cl-]' });
            if (lower.includes('coffee')) chemicals.push({ name: 'caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' });
            if (lower.includes('apple')) chemicals.push({ name: 'fructose', smiles: 'C(C(C(C(C(C=O)O)O)O)O)O' });
            if (chemicals.length === 0) chemicals.push({ name: 'water', smiles: 'O' });

            const json = { object: text || 'test object', chemicals };
            return { choices: [{ message: { content: JSON.stringify(json) } }] };
          })
        }
      }
    }))
  };
});

// Set test environment variables (do not force OPENAI_API_KEY; allow skipping heavy tests)
process.env.NODE_ENV = "test";
if (process.env.OPENAI_API_KEY) {
  delete process.env.OPENAI_API_KEY;
}

// Polyfill setImmediate for test environment
if (!global.setImmediate) {
  global.setImmediate = (callback, ...args) => setTimeout(callback, 0, ...args);
}

// Setup process cleanup hooks
process.on('exit', () => cleanup.forceCleanup());
process.on('SIGINT', async () => {
  await cleanup.cleanup();
  process.exit(0);
});
process.on('SIGTERM', async () => {
  await cleanup.cleanup();
  process.exit(0);
});

// Override global timer functions to track them
const originalSetTimeout = global.setTimeout;
const originalSetInterval = global.setInterval;

global.setTimeout = (...args) => {
  const timer = originalSetTimeout(...args);
  cleanup.registerTimer(timer);
  return timer;
};

global.setInterval = (...args) => {
  const interval = originalSetInterval(...args);
  cleanup.registerInterval(interval);
  return interval;
};
