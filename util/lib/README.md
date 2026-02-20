# Inline Debugger Library

A comprehensive debugging library for JavaScript/Node.js that provides inline inspection, breakpoint-like functionality, performance profiling, and advanced debugging capabilities.

## Features

- 🔍 **Inline Variable Inspection** - Inspect variables with detailed type and size information
- 🛑 **Breakpoint-like Functionality** - Conditional stopping with context information
- 📊 **Performance Profiling** - Time functions and monitor memory usage
- 👀 **Variable Watching** - Monitor variable changes in real-time
- 📍 **Function Tracing** - Track function calls with data and stack traces
- 🧠 **Memory Monitoring** - Monitor memory usage and garbage collection
- 📝 **Comprehensive Logging** - Multiple log levels with file output support
- ⚡ **Performance Monitoring** - Built-in performance measurement tools

## Installation

```bash
# Copy the library to your project
cp -r util/lib/debug.js ./lib/
```

## Quick Start

```javascript
const debug = require('./lib/debug');

// Basic inspection
const data = { name: 'John', age: 30 };
debug.inspect(data, 'User Data');

// Breakpoint with condition
debug.breakpoint('checkpoint', data.age > 25, 'User is over 25');

// Function profiling
const profiledFunction = debug.profile('myFunction', myFunction);

// Timing operations
debug.time('operation');
// ... do work ...
debug.timeEnd('operation');

// Memory monitoring
debug.memory();
```

## API Reference

### Core Methods

#### `debug.inspect(variable, label, options)`
Inspect a variable with detailed information.

```javascript
debug.inspect(myObject, 'My Object', { 
  depth: 3, 
  showType: true, 
  showSize: true 
});
```

#### `debug.breakpoint(name, condition, message)`
Create a conditional breakpoint.

```javascript
debug.breakpoint('user-check', user.isAdmin, 'Admin user detected');
debug.breakpoint('counter', () => counter > 10, 'Counter exceeded limit');
```

#### `debug.trace(functionName, data, level)`
Trace function execution.

```javascript
function myFunction(data) {
  debug.trace('myFunction', data);
  // ... function logic
}
```

### Timing and Performance

#### `debug.time(label)` / `debug.timeEnd(label)`
Time operations.

```javascript
debug.time('data-processing');
// ... do work ...
const duration = debug.timeEnd('data-processing');
```

#### `debug.profile(functionName, fn)`
Profile function performance.

```javascript
const profiledFn = debug.profile('expensiveOperation', expensiveOperation);
```

#### `debug.perf(label, fn)`
Performance monitoring wrapper.

```javascript
const monitoredFn = debug.perf('myFunction', myFunction);
```

### Advanced Features

#### `debug.watch(variableName, getter, options)`
Watch variable changes.

```javascript
let myVar = 0;
debug.watch('myVar', () => myVar, { interval: 1000 });
```

#### `debug.dump(variable, filename)`
Dump variable to file or console.

```javascript
debug.dump(complexObject, './debug-dump.json');
```

#### `debug.memory()`
Monitor memory usage.

```javascript
debug.memory(); // Logs current memory usage
```

### Logging

#### Log Levels
```javascript
debug.info('Info message', data);
debug.warn('Warning message', data);
debug.error('Error message', data);
debug.success('Success message', data);
debug.debug('Debug message', data);
```

### Utility Methods

#### `debug.clear()`
Clear all debug state.

```javascript
debug.clear(); // Clears timers, breakpoints, traces, watchers
```

#### `debug.getStats()`
Get debug statistics.

```javascript
const stats = debug.getStats();
console.log(stats);
```

## Configuration

### Default Instance
```javascript
const debug = require('./lib/debug');
```

### Custom Instance
```javascript
const { DebugLib } = require('./lib/debug');
const customDebug = new DebugLib({
  enabled: true,
  logLevel: 'debug',
  outputFile: './debug.log',
  depth: 5,
  colors: true,
  maxArrayLength: 100
});
```

### Environment Variables
- `NODE_ENV=production` - Disables debugging
- `DEBUG_LOG_FILE=./debug.log` - Sets log file path

## Examples

### Basic Usage
```javascript
const debug = require('./lib/debug');

// Inspect variables
const user = { name: 'John', age: 30 };
debug.inspect(user, 'User Object');

// Time operations
debug.time('processing');
// ... do work ...
debug.timeEnd('processing');

// Conditional breakpoints
debug.breakpoint('user-check', user.age > 18, 'User is adult');
```

### Advanced Usage
```javascript
// Function profiling
const expensiveFunction = debug.profile('expensive', (data) => {
  // ... expensive operation
});

// Variable watching
let counter = 0;
debug.watch('counter', () => counter, { interval: 1000 });

// Performance monitoring
const monitoredFunction = debug.perf('myFunction', myFunction);

// Memory monitoring
setInterval(() => debug.memory(), 5000);
```

## Running Examples

```bash
# Basic usage examples
node examples/basic-usage.js

# Advanced usage examples
node examples/advanced-usage.js
```

## License

MIT