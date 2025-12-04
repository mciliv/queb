/**
 * Common test helper functions
 */

const fs = require('fs');
const path = require('path');

/**
 * Sleep for specified milliseconds
 */
const sleep = (ms) => new Promise(resolve => setTimeout(resolve, ms));

/**
 * Retry a function until it succeeds or timeout
 */
const retry = async (fn, options = {}) => {
  const {
    retries = 3,
    delay = 1000,
    timeout = 10000,
    onRetry = () => {}
  } = options;
  
  const startTime = Date.now();
  let lastError;
  
  for (let i = 0; i < retries; i++) {
    try {
      return await fn();
    } catch (error) {
      lastError = error;
      
      if (Date.now() - startTime > timeout) {
        throw new Error(`Timeout after ${timeout}ms: ${error.message}`);
      }
      
      if (i < retries - 1) {
        onRetry(i + 1, error);
        await sleep(delay);
      }
    }
  }
  
  throw lastError;
};

/**
 * Create a temporary test file
 */
const createTempFile = (content, extension = 'txt') => {
  const tempDir = path.join(__dirname, '../__temp__');
  if (!fs.existsSync(tempDir)) {
    fs.mkdirSync(tempDir, { recursive: true });
  }
  
  const filename = `test-${Date.now()}-${Math.random().toString(36).substr(2)}.${extension}`;
  const filepath = path.join(tempDir, filename);
  
  fs.writeFileSync(filepath, content);
  
  return {
    path: filepath,
    cleanup: () => {
      if (fs.existsSync(filepath)) {
        fs.unlinkSync(filepath);
      }
    }
  };
};

/**
 * Clean up all temporary test files
 */
const cleanupTempFiles = () => {
  const tempDir = path.join(__dirname, '../__temp__');
  if (fs.existsSync(tempDir)) {
    fs.rmSync(tempDir, { recursive: true, force: true });
  }
};

/**
 * Mock HTTP response
 */
const mockResponse = (data, options = {}) => {
  return {
    ok: options.ok !== false,
    status: options.status || 200,
    statusText: options.statusText || 'OK',
    headers: new Map(Object.entries(options.headers || {})),
    json: async () => data,
    text: async () => JSON.stringify(data),
    blob: async () => new Blob([JSON.stringify(data)], { type: 'application/json' }),
    clone: function() { return { ...this }; }
  };
};

/**
 * Create a mock API client
 */
const createMockApiClient = () => {
  const responses = new Map();
  const requests = [];
  
  return {
    // Set up a mock response
    mockEndpoint: (method, url, response) => {
      const key = `${method.toUpperCase()} ${url}`;
      responses.set(key, response);
    },
    
    // Get all requests made
    getRequests: () => requests,
    
    // Clear all mocks and requests
    reset: () => {
      responses.clear();
      requests.length = 0;
    },
    
    // Mock fetch function
    fetch: async (url, options = {}) => {
      const method = options.method || 'GET';
      const key = `${method} ${url}`;
      
      requests.push({ url, ...options, timestamp: Date.now() });
      
      if (responses.has(key)) {
        const response = responses.get(key);
        if (typeof response === 'function') {
          return response(url, options);
        }
        return mockResponse(response);
      }
      
      throw new Error(`No mock response for ${key}`);
    }
  };
};

/**
 * Wait for a condition to be true
 */
const waitForCondition = async (condition, message = 'Condition not met', timeout = 5000) => {
  const interval = 100;
  const startTime = Date.now();
  
  while (Date.now() - startTime < timeout) {
    if (await condition()) {
      return;
    }
    await sleep(interval);
  }
  
  throw new Error(`${message} after ${timeout}ms`);
};

/**
 * Create test data generators
 */
const generators = {
  // Generate random ID
  id: () => Math.random().toString(36).substr(2, 9),
  
  // Generate random email
  email: () => `test-${Date.now()}@example.com`,
  
  // Generate random SMILES
  smiles: () => {
    const smiles = ['CCO', 'CC(=O)O', 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O', 'O', 'C1=CC=CC=C1'];
    return smiles[Math.floor(Math.random() * smiles.length)];
  },
  
  // Generate random molecule name
  moleculeName: () => {
    const names = ['Ethanol', 'Acetic acid', 'Ibuprofen', 'Water', 'Benzene'];
    return names[Math.floor(Math.random() * names.length)];
  },
  
  // Generate molecular data
  molecule: () => ({
    name: generators.moleculeName(),
    smiles: generators.smiles(),
    id: generators.id()
  }),
  
  // Generate timestamp
  timestamp: () => new Date().toISOString()
};

/**
 * Performance measurement helper
 */
class PerformanceMeasure {
  constructor(name) {
    this.name = name;
    this.startTime = Date.now();
    this.marks = [];
  }
  
  mark(label) {
    const elapsed = Date.now() - this.startTime;
    this.marks.push({ label, elapsed });
    return elapsed;
  }
  
  end() {
    const totalTime = Date.now() - this.startTime;
    return {
      name: this.name,
      totalTime,
      marks: this.marks
    };
  }
}

/**
 * Test data validator
 */
const validators = {
  isValidSmiles: (smiles) => {
    return typeof smiles === 'string' && /^[A-Za-z0-9@+\-\[\]()=#$\\.\\/:]+$/.test(smiles);
  },
  
  isValidMolecule: (molecule) => {
    return molecule &&
      typeof molecule === 'object' &&
      typeof molecule.name === 'string' &&
      validators.isValidSmiles(molecule.smiles);
  },
  
  isValidUrl: (url) => {
    try {
      new URL(url);
      return true;
    } catch {
      return false;
    }
  }
};

/**
 * Mock console for capturing logs
 */
class MockConsole {
  constructor() {
    this.logs = [];
    this.originalConsole = { ...console };
  }
  
  start() {
    ['log', 'info', 'warn', 'error', 'debug'].forEach(method => {
      console[method] = (...args) => {
        this.logs.push({ method, args, timestamp: Date.now() });
      };
    });
  }
  
  stop() {
    Object.assign(console, this.originalConsole);
  }
  
  getLogs(method) {
    return method 
      ? this.logs.filter(log => log.method === method)
      : this.logs;
  }
  
  clear() {
    this.logs = [];
  }
}

module.exports = {
  sleep,
  retry,
  createTempFile,
  cleanupTempFiles,
  mockResponse,
  createMockApiClient,
  waitForCondition,
  generators,
  PerformanceMeasure,
  validators,
  MockConsole
};
