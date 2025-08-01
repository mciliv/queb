// Production-safe logger - replaces direct console.log usage
const debug = const IS_PRODUCTION = process.env.NODE_ENV === 'production';
const LOG_LEVEL = process.env.LOG_LEVEL || (IS_PRODUCTION ? 'error' : 'debug');

const LOG_LEVELS = {
  debug: 0,
  info: 1, 
  warn: 2,
  error: 3
};

class Logger {
  constructor(module = 'APP') {
    this.module = module;
  }

  debug(...args) {
    if (this.shouldLog('debug')) {
      
    }
  }

  info(...args) {
    if (this.shouldLog('info')) {
      console.info(`[${this.module}]`, ...args);
    }
  }

  warn(...args) {
    if (this.shouldLog('warn')) {
      
    }
  }

  error(...args) {
    if (this.shouldLog('error')) {
      console.error(`[${this.module}]`, ...args);
    }
  }

  // Time tracking (debug only)
  time(label) {
    debug.time(`${this.module}.${label}`);
  }

  timeEnd(label) {
    debug.timeEnd(`${this.module}.${label}`);
  }

  // Object inspection (debug only)
  inspect(obj, label) {
    
  }

  shouldLog(level) {
    return LOG_LEVELS[level] >= LOG_LEVELS[LOG_LEVEL];
  }
}

// Singleton instances for common modules
const loggers = {
  server: new Logger('SERVER'),
  api: new Logger('API'),
  db: new Logger('DATABASE'),
  auth: new Logger('AUTH'),
  molecules: new Logger('MOLECULES'),
  payment: new Logger('PAYMENT')
};

module.exports = {
  Logger,
  ...loggers
}; 