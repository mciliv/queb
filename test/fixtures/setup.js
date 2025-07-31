// test/fixtures/setup.js - Jest setup file for mocking dependencies

const cleanup = require('./cleanup-registry');

// Mock OpenAI before any modules are imported
jest.mock("openai", () => {
  return {
    OpenAI: jest.fn().mockImplementation(() => ({
      responses: {
        parse: jest.fn().mockResolvedValue({
          output_parsed: {
            object: "test object",
            smiles: ["O", "CCO"],
          },
        }),
      },
    })),
  };
});

// Set test environment variables
process.env.NODE_ENV = "test";
process.env.OPENAI_API_KEY = "test-key";

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
