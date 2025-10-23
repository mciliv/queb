#!/usr/bin/env node

/**
 * Inline Debugger Utility
 * The closest thing to an inline debugger for JavaScript/Node.js
 * 
 * Usage:
 *   const debug = require('./debug');
 *   debug.inspect(variable);
 *   debug.breakpoint('checkpoint1', condition);
 *   debug.trace('function_name', data);
 *   debug.time('operation');
 */

const util = require('util');
const fs = require('fs');
const path = require('path');

class InlineDebugger {
  constructor(options = {}) {
    this.enabled = options.enabled !== false;
    this.logLevel = options.logLevel || 'info';
    this.outputFile = options.outputFile || null;
    this.colors = options.colors !== false;
    this.depth = options.depth || 3;
    this.maxArrayLength = options.maxArrayLength || 100;
    this.timers = new Map();
    this.breakpoints = new Map();
    this.traces = [];
    
    // Color codes for console output
    this.colors = {
      reset: '\x1b[0m',
      bright: '\x1b[1m',
      dim: '\x1b[2m',
      red: '\x1b[31m',
      green: '\x1b[32m',
      yellow: '\x1b[33m',
      blue: '\x1b[34m',
      magenta: '\x1b[35m',
      cyan: '\x1b[36m',
      white: '\x1b[37m',
      bgRed: '\x1b[41m',
      bgGreen: '\x1b[42m',
      bgYellow: '\x1b[43m',
      bgBlue: '\x1b[44m'
    };
  }

  _colorize(text, color) {
    if (!this.colors) return text;
    return `${this.colors[color]}${text}${this.colors.reset}`;
  }

  _log(level, message, data = null) {
    if (!this.enabled) return;
    
    const timestamp = new Date().toISOString();
    const levelColor = {
      'info': 'cyan',
      'warn': 'yellow', 
      'error': 'red',
      'success': 'green',
      'debug': 'dim'
    }[level] || 'white';
    
    const prefix = this._colorize(`[${timestamp}] [${level.toUpperCase()}]`, levelColor);
    const output = data ? `${prefix} ${message}\n${this._inspect(data)}` : `${prefix} ${message}`;
    
    console.log(output);
    
    if (this.outputFile) {
      this._writeToFile(level, message, data);
    }
  }

  _writeToFile(level, message, data) {
    try {
      const logEntry = {
        timestamp: new Date().toISOString(),
        level,
        message,
        data: data ? this._serialize(data) : null
      };
      
      const logLine = JSON.stringify(logEntry) + '\n';
      fs.appendFileSync(this.outputFile, logLine);
    } catch (error) {
      console.error('Failed to write to debug log file:', error.message);
    }
  }

  _serialize(obj) {
    const seen = new WeakSet();
    return JSON.parse(JSON.stringify(obj, (key, value) => {
      if (value instanceof Error) {
        return { message: value.message, stack: value.stack, name: value.name };
      }
      if (typeof value === 'object' && value !== null) {
        if (seen.has(value)) return '[Circular Reference]';
        seen.add(value);
      }
      return value;
    }));
  }

  _inspect(obj, customDepth = null) {
    const depth = customDepth || this.depth;
    const inspectOptions = {
      depth,
      colors: this.colors,
      maxArrayLength: this.maxArrayLength,
      showHidden: true,
      compact: false
    };
    
    return util.inspect(obj, inspectOptions);
  }

  // Main debugging methods
  inspect(variable, label = 'Variable', options = {}) {
    if (!this.enabled) return variable;
    
    const { depth, showType = true, showSize = true } = options;
    
    let output = this._colorize(`\n🔍 ${label}:`, 'bright');
    
    if (showType) {
      const type = Array.isArray(variable) ? 'Array' : typeof variable;
      const size = Array.isArray(variable) ? variable.length : 
                   typeof variable === 'object' && variable !== null ? Object.keys(variable).length : 'N/A';
      
      output += `\n   Type: ${this._colorize(type, 'blue')}`;
      if (showSize && size !== 'N/A') {
        output += ` | Size: ${this._colorize(size, 'green')}`;
      }
    }
    
    output += `\n${this._inspect(variable, depth)}`;
    
    console.log(output);
    return variable;
  }

  breakpoint(name, condition = true, message = null) {
    if (!this.enabled) return;
    
    const shouldBreak = typeof condition === 'function' ? condition() : condition;
    
    if (shouldBreak) {
      const breakMsg = message || `Breakpoint: ${name}`;
      this._log('warn', `🛑 ${breakMsg}`);
      
      // In a real debugger, this would pause execution
      // Here we provide detailed context
      console.log(this._colorize('\n📊 Breakpoint Context:', 'bright'));
      console.log('   Call Stack:', new Error().stack.split('\n').slice(1, 6).join('\n'));
      console.log('   Memory Usage:', process.memoryUsage());
      
      this.breakpoints.set(name, {
        hit: true,
        timestamp: new Date().toISOString(),
        message: breakMsg
      });
    }
  }

  trace(functionName, data = null, level = 'info') {
    if (!this.enabled) return;
    
    const traceEntry = {
      function: functionName,
      timestamp: new Date().toISOString(),
      data: data ? this._serialize(data) : null,
      memory: process.memoryUsage(),
      stack: new Error().stack.split('\n').slice(1, 4)
    };
    
    this.traces.push(traceEntry);
    this._log(level, `📍 Trace: ${functionName}`, data);
  }

  time(label) {
    if (!this.enabled) return;
    
    if (this.timers.has(label)) {
      const startTime = this.timers.get(label);
      const duration = Date.now() - startTime;
      this._log('info', `⏱️  Timer '${label}': ${duration}ms`);
      this.timers.delete(label);
      return duration;
    } else {
      this.timers.set(label, Date.now());
      this._log('debug', `⏱️  Timer '${label}' started`);
    }
  }

  timeEnd(label) {
    return this.time(label);
  }

  // Advanced debugging methods
  dump(variable, filename = null) {
    if (!this.enabled) return;
    
    const dumpData = {
      timestamp: new Date().toISOString(),
      variable: this._serialize(variable),
      memory: process.memoryUsage(),
      stack: new Error().stack
    };
    
    if (filename) {
      const filepath = path.resolve(filename);
      fs.writeFileSync(filepath, JSON.stringify(dumpData, null, 2));
      this._log('info', `💾 Variable dumped to: ${filepath}`);
    } else {
      console.log(this._colorize('\n💾 Variable Dump:', 'bright'));
      console.log(JSON.stringify(dumpData, null, 2));
    }
  }

  profile(functionName, fn) {
    if (!this.enabled) return fn;
    
    return (...args) => {
      const startTime = process.hrtime.bigint();
      const startMemory = process.memoryUsage();
      
      try {
        const result = fn.apply(this, args);
        const endTime = process.hrtime.bigint();
        const endMemory = process.memoryUsage();
        
        const duration = Number(endTime - startTime) / 1000000; // Convert to milliseconds
        const memoryDelta = {
          rss: endMemory.rss - startMemory.rss,
          heapUsed: endMemory.heapUsed - startMemory.heapUsed,
          heapTotal: endMemory.heapTotal - startMemory.heapTotal
        };
        
        this._log('info', `📊 Profile: ${functionName}`, {
          duration: `${duration.toFixed(2)}ms`,
          memoryDelta: `${Math.round(memoryDelta.heapUsed / 1024)}KB`,
          args: args.length,
          result: typeof result
        });
        
        return result;
      } catch (error) {
        this._log('error', `❌ Profile Error: ${functionName}`, error.message);
        throw error;
      }
    };
  }

  // Utility methods
  clear() {
    this.timers.clear();
    this.breakpoints.clear();
    this.traces = [];
    this._log('info', '🧹 Debug state cleared');
  }

  getStats() {
    return {
      timers: this.timers.size,
      breakpoints: this.breakpoints.size,
      traces: this.traces.length,
      memory: process.memoryUsage()
    };
  }

  // Console methods for different log levels
  info(message, data = null) {
    this._log('info', message, data);
  }

  warn(message, data = null) {
    this._log('warn', message, data);
  }

  error(message, data = null) {
    this._log('error', message, data);
  }

  success(message, data = null) {
    this._log('success', message, data);
  }

  debug(message, data = null) {
    this._log('debug', message, data);
  }
}

// Create default instance
const debug = new InlineDebugger({
  enabled: process.env.NODE_ENV !== 'production',
  outputFile: process.env.DEBUG_LOG_FILE || null
});

// Export both the class and default instance
module.exports = debug;
module.exports.InlineDebugger = InlineDebugger;

// If running directly, provide CLI interface
if (require.main === module) {
  const args = process.argv.slice(2);
  
  if (args.length === 0) {
    console.log(`
🔍 Inline Debugger Utility

Usage:
  node debug.js inspect <file> <line> <variable>
  node debug.js profile <function>
  node debug.js stats
  node debug.js clear

Examples:
  node debug.js inspect ./app.js 42 myVariable
  node debug.js profile myFunction
  node debug.js stats
    `);
    process.exit(0);
  }
  
  const command = args[0];
  
  switch (command) {
    case 'stats':
      console.log('Debug Statistics:', debug.getStats());
      break;
    case 'clear':
      debug.clear();
      break;
    default:
      console.log('Unknown command:', command);
  }
}