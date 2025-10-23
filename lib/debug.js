/**
 * Inline Debugger Library
 * A comprehensive debugging library for JavaScript/Node.js
 * 
 * Features:
 * - Inline variable inspection
 * - Breakpoint-like functionality
 * - Performance profiling
 * - Memory monitoring
 * - Function tracing
 * - Conditional debugging
 */

const util = require('util');
const fs = require('fs');
const path = require('path');

class DebugLib {
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
    this.watchers = new Map();
    
    // Color codes for console output
    this.colorCodes = {
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
    return `${this.colorCodes[color]}${text}${this.colorCodes.reset}`;
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

  // Core debugging methods
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
        
        const duration = Number(endTime - startTime) / 1000000;
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

  // Variable watching
  watch(variableName, getter, options = {}) {
    if (!this.enabled) return;
    
    const { interval = 1000, onChange = true } = options;
    let lastValue = null;
    
    const watcher = setInterval(() => {
      try {
        const currentValue = typeof getter === 'function' ? getter() : getter;
        
        if (onChange && JSON.stringify(currentValue) !== JSON.stringify(lastValue)) {
          this._log('info', `👀 Watch: ${variableName} changed`, {
            from: lastValue,
            to: currentValue
          });
          lastValue = currentValue;
        }
      } catch (error) {
        this._log('error', `❌ Watch Error: ${variableName}`, error.message);
      }
    }, interval);
    
    this.watchers.set(variableName, watcher);
    this._log('debug', `👀 Started watching: ${variableName}`);
  }

  unwatch(variableName) {
    const watcher = this.watchers.get(variableName);
    if (watcher) {
      clearInterval(watcher);
      this.watchers.delete(variableName);
      this._log('debug', `👀 Stopped watching: ${variableName}`);
    }
  }

  // Memory monitoring
  memory() {
    if (!this.enabled) return;
    
    const mem = process.memoryUsage();
    const format = (bytes) => `${Math.round(bytes / 1024 / 1024 * 100) / 100} MB`;
    
    this._log('info', '🧠 Memory Usage:', {
      rss: format(mem.rss),
      heapTotal: format(mem.heapTotal),
      heapUsed: format(mem.heapUsed),
      external: format(mem.external),
      arrayBuffers: format(mem.arrayBuffers)
    });
  }

  // Performance monitoring
  perf(label, fn) {
    if (!this.enabled) return fn;
    
    return (...args) => {
      const start = process.hrtime.bigint();
      const startMem = process.memoryUsage();
      
      try {
        const result = fn.apply(this, args);
        const end = process.hrtime.bigint();
        const endMem = process.memoryUsage();
        
        const duration = Number(end - start) / 1000000;
        const memoryUsed = endMem.heapUsed - startMem.heapUsed;
        
        this._log('info', `⚡ Performance: ${label}`, {
          duration: `${duration.toFixed(2)}ms`,
          memoryUsed: `${Math.round(memoryUsed / 1024)}KB`,
          args: args.length
        });
        
        return result;
      } catch (error) {
        this._log('error', `❌ Performance Error: ${label}`, error.message);
        throw error;
      }
    };
  }

  // Utility methods
  clear() {
    this.timers.clear();
    this.breakpoints.clear();
    this.traces = [];
    this.watchers.forEach(watcher => clearInterval(watcher));
    this.watchers.clear();
    this._log('info', '🧹 Debug state cleared');
  }

  getStats() {
    return {
      timers: this.timers.size,
      breakpoints: this.breakpoints.size,
      traces: this.traces.length,
      watchers: this.watchers.size,
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
const debug = new DebugLib({
  enabled: process.env.NODE_ENV !== 'production',
  outputFile: process.env.DEBUG_LOG_FILE || null
});

// Export both the class and default instance
module.exports = debug;
module.exports.DebugLib = DebugLib;