#!/usr/bin/env node

/**
 * Advanced Usage Examples for Inline Debugger
 */

const debug = require('../debug');

console.log('🔍 Inline Debugger - Advanced Usage Examples\n');

// 1. Function profiling
function fibonacci(n) {
  if (n <= 1) return n;
  return fibonacci(n - 1) + fibonacci(n - 2);
}

const profiledFibonacci = debug.profile('fibonacci', fibonacci);
const result = profiledFibonacci(10);
debug.inspect(result, 'Fibonacci Result');

// 2. Performance monitoring
const slowOperation = debug.perf('slow-operation', (data) => {
  // Simulate slow operation
  let sum = 0;
  for (let i = 0; i < 1000000; i++) {
    sum += Math.random();
  }
  return sum;
});

slowOperation('test-data');

// 3. Variable watching
let watchedVariable = 0;
debug.watch('watchedVariable', () => watchedVariable, { interval: 500 });

// Change the variable
setTimeout(() => {
  watchedVariable = 42;
}, 1000);

setTimeout(() => {
  watchedVariable = 100;
}, 2000);

// 4. Data dumping
const complexData = {
  users: [
    { id: 1, name: 'Alice', scores: [85, 92, 78] },
    { id: 2, name: 'Bob', scores: [91, 88, 95] }
  ],
  metadata: {
    total: 2,
    average: 88.5,
    timestamp: new Date()
  }
};

debug.dump(complexData, './debug-dump.json');

// 5. Custom debug instance with options
const { DebugLib } = require('../debug');
const customDebug = new DebugLib({
  enabled: true,
  logLevel: 'debug',
  outputFile: './custom-debug.log',
  depth: 5,
  colors: true
});

customDebug.info('Custom debug instance message');
customDebug.inspect({ custom: 'data' }, 'Custom Debug Data');

// 6. Error handling and debugging
function riskyOperation(data) {
  try {
    if (!data) {
      throw new Error('Data is required');
    }
    
    debug.trace('riskyOperation', { data });
    
    // Simulate some processing
    const result = data.map(item => item * 2);
    
    debug.success('Risky operation completed successfully');
    return result;
  } catch (error) {
    debug.error('Risky operation failed', error.message);
    throw error;
  }
}

// Test with valid data
try {
  const validData = [1, 2, 3, 4, 5];
  const result = riskyOperation(validData);
  debug.inspect(result, 'Valid Operation Result');
} catch (error) {
  debug.error('Operation failed', error.message);
}

// Test with invalid data
try {
  riskyOperation(null);
} catch (error) {
  debug.error('Expected error caught', error.message);
}

// 7. Memory monitoring over time
setInterval(() => {
  debug.memory();
}, 2000);

// 8. Debug statistics
setTimeout(() => {
  const stats = debug.getStats();
  debug.inspect(stats, 'Debug Statistics');
  
  // Clean up
  debug.unwatch('watchedVariable');
  debug.clear();
  
  console.log('\n✅ Advanced usage examples completed!');
  process.exit(0);
}, 5000);