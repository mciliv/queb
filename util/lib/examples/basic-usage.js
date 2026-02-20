#!/usr/bin/env node

/**
 * Basic Usage Examples for Inline Debugger
 */

const debug = require('../debug');

console.log('🔍 Inline Debugger - Basic Usage Examples\n');

// 1. Basic variable inspection
const user = {
  name: 'John Doe',
  age: 30,
  email: 'john@example.com',
  preferences: {
    theme: 'dark',
    notifications: true
  }
};

debug.inspect(user, 'User Object');

// 2. Array inspection
const numbers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
debug.inspect(numbers, 'Numbers Array');

// 3. Function tracing
function calculateSum(a, b) {
  debug.trace('calculateSum', { a, b });
  return a + b;
}

const result = calculateSum(5, 3);
debug.inspect(result, 'Sum Result');

// 4. Timing operations
debug.time('data-processing');

// Simulate some work
setTimeout(() => {
  debug.timeEnd('data-processing');
}, 100);

// 5. Breakpoint with condition
let counter = 0;
setInterval(() => {
  counter++;
  debug.breakpoint('counter-check', counter > 3, `Counter reached ${counter}`);
  
  if (counter >= 5) {
    process.exit(0);
  }
}, 200);

// 6. Memory monitoring
debug.memory();

// 7. Different log levels
debug.info('This is an info message');
debug.warn('This is a warning message');
debug.success('This is a success message');
debug.error('This is an error message');
debug.debug('This is a debug message');

console.log('\n✅ Basic usage examples completed!');