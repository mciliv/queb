#!/usr/bin/env node

// Get the entry point from command line arguments
const entryPoint = process.argv[2];

if (!entryPoint) {
  console.error('Error: Entry point required as argument');
  process.exit(1);
}

// Polyfill browser globals for server environment
global.window = global.window || {
  location: { hostname: 'localhost', search: '', href: 'http://localhost:8080' },
  navigator: { userAgent: 'Node.js Server' },
  localStorage: { 
    getItem: () => null, 
    setItem: () => {}, 
    removeItem: () => {} 
  },
  addEventListener: () => {},
  removeEventListener: () => {}
};

global.document = global.document || {
  addEventListener: () => {},
  removeEventListener: () => {}
};

global.localStorage = global.window.localStorage;

// Now require and start the actual server
try {
  require(entryPoint);
} catch (error) {
  console.error('Failed to start server:', error.message);
  process.exit(1);
}



