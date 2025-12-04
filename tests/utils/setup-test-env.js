/**
 * Test environment setup
 * This runs before each test file
 */

// Polyfills for Node.js environment
global.TextEncoder = TextEncoder;
global.TextDecoder = TextDecoder;

// Mock console methods to reduce noise during tests
const originalConsole = { ...console };
global.console = {
  ...console,
  log: jest.fn(),
  info: jest.fn(),
  warn: jest.fn(),
  // Keep error and debug for important messages
  error: originalConsole.error,
  debug: process.env.DEBUG ? originalConsole.debug : jest.fn()
};

// Add custom Jest matchers
expect.extend({
  // Check if a value is a valid SMILES string
  toBeValidSmiles(received) {
    const smilesRegex = /^[A-Za-z0-9@+\-\[\]()=#$\\.\\/:]+$/;
    const pass = typeof received === 'string' && smilesRegex.test(received);
    
    if (pass) {
      return {
        message: () => `expected ${received} not to be a valid SMILES string`,
        pass: true
      };
    } else {
      return {
        message: () => `expected ${received} to be a valid SMILES string`,
        pass: false
      };
    }
  },

  // Check if an object has expected molecular structure
  toBeMolecularData(received) {
    const pass = 
      received &&
      typeof received === 'object' &&
      'name' in received &&
      'smiles' in received;
    
    if (pass) {
      return {
        message: () => `expected ${JSON.stringify(received)} not to be molecular data`,
        pass: true
      };
    } else {
      return {
        message: () => `expected ${JSON.stringify(received)} to be molecular data with 'name' and 'smiles' properties`,
        pass: false
      };
    }
  },

  // Check if a response is successful
  toBeSuccessResponse(received) {
    const pass = 
      received &&
      typeof received === 'object' &&
      (received.success === true || received.status === 'success' || received.ok === true);
    
    if (pass) {
      return {
        message: () => `expected response not to be successful`,
        pass: true
      };
    } else {
      return {
        message: () => `expected response to be successful (have success: true, status: 'success', or ok: true)`,
        pass: false
      };
    }
  }
});

// Common test utilities
global.testUtils = {
  // Wait for a condition to be true
  waitFor: async (condition, timeout = 5000, interval = 100) => {
    const startTime = Date.now();
    while (Date.now() - startTime < timeout) {
      if (await condition()) {
        return true;
      }
      await new Promise(resolve => setTimeout(resolve, interval));
    }
    throw new Error(`Timeout waiting for condition after ${timeout}ms`);
  },

  // Create a deferred promise
  createDeferred: () => {
    let resolve, reject;
    const promise = new Promise((res, rej) => {
      resolve = res;
      reject = rej;
    });
    return { promise, resolve, reject };
  },

  // Mock timer that can be controlled
  createMockTimer: () => {
    const callbacks = [];
    return {
      setTimeout: (cb, delay) => {
        const id = callbacks.length;
        callbacks.push({ cb, delay, id });
        return id;
      },
      clearTimeout: (id) => {
        const index = callbacks.findIndex(c => c.id === id);
        if (index !== -1) callbacks.splice(index, 1);
      },
      tick: (ms) => {
        const toRun = callbacks.filter(c => c.delay <= ms);
        callbacks.splice(0, callbacks.length, ...callbacks.filter(c => c.delay > ms));
        toRun.forEach(c => c.cb());
      },
      clear: () => callbacks.splice(0, callbacks.length)
    };
  }
};

// Clean up after each test
afterEach(() => {
  // Clear all mocks
  jest.clearAllMocks();
  
  // Clear any test timeouts
  jest.clearAllTimers();
});
