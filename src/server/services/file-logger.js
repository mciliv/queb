const fs = require('fs');
const path = require('path');

class FileLogger {
  constructor() {
    // Use the project root for logs, not the current working directory
    const projectRoot = this.findProjectRoot();
    this.logsDir = process.env.LOGS_DIR || path.join(projectRoot, 'logs');
    this.minimumLevel = (process.env.LOG_LEVEL || 'info').toLowerCase();
    this.levelOrder = { error: 0, warn: 1, success: 2, info: 2, log: 2, debug: 3 };
    this.ensureLogsDirectory();
    this.clearOldLogs();
  }

  findProjectRoot() {
    // Start from the current file location and walk up to find project root
    let currentDir = __dirname;
    
    // Walk up the directory tree to find package.json (project root indicator)
    while (currentDir !== path.dirname(currentDir)) {
      if (fs.existsSync(path.join(currentDir, 'package.json'))) {
        return currentDir;
      }
      currentDir = path.dirname(currentDir);
    }
    
    // If we're in a src/ directory, go up one level to find the project root
    if (currentDir.endsWith('src')) {
      const parentDir = path.dirname(currentDir);
      if (fs.existsSync(path.join(parentDir, 'package.json'))) {
        return parentDir;
      }
    }
    
    // Fallback to process.cwd() if we can't find package.json
    return process.cwd();
  }

  ensureLogsDirectory() {
    try {
      if (!fs.existsSync(this.logsDir)) {
        fs.mkdirSync(this.logsDir, { recursive: true });
      }
    } catch (error) {
      console.warn('⚠️ Could not create logs directory, falling back to console logging:', error.message);
      this.logsDir = null;
    }
  }

  clearOldLogs() {
    if (!this.logsDir) return;
    
    try {
      const files = fs.readdirSync(this.logsDir);
      let clearedCount = 0;

      for (const file of files) {
        // Only delete .log files, preserve screenshots and other files
        if (file.endsWith('.log') && !file.includes('Screenshot')) {
          const filePath = path.join(this.logsDir, file);
          fs.unlinkSync(filePath);
          clearedCount++;
        }
      }

      if (clearedCount > 0) {
        console.log(`🧹 Cleared ${clearedCount} old log file(s) from previous session`);
        // Also write to the new log file
        this.writeToFile(this.getLogFile(), this.formatFileLine('info', `🧹 Cleared ${clearedCount} old log file(s) from previous session`));
      }
    } catch (error) {
      console.warn('⚠️ Failed to clear old logs:', error.message);
    }
  }

  getLogFile(type = 'server') {
    if (!this.logsDir) return null;
    const today = new Date().toISOString().split('T')[0];
    return path.join(this.logsDir, `${type}-${today}.log`);
  }

  safeSerializeMeta(meta) {
    if (!meta || typeof meta !== 'object') return undefined;
    const seen = new WeakSet();
    return JSON.parse(JSON.stringify(meta, (key, value) => {
      if (value instanceof Error) {
        return { message: value.message, stack: value.stack, name: value.name };
      }
      if (typeof value === 'object' && value !== null) {
        if (seen.has(value)) return '[Circular]';
        seen.add(value);
      }
      return value;
    }));
  }

  shouldLog(level) {
    const current = this.levelOrder[level] ?? 2;
    const min = this.levelOrder[this.minimumLevel] ?? 2;
    return current <= min;
  }

  formatConsole(level, message, meta) {
    const timestamp = new Date().toISOString();
    const levelLabel = level.toUpperCase();
    const color = level === 'error' ? '\x1b[31m' : level === 'warn' ? '\x1b[33m' : level === 'debug' ? '\x1b[90m' : level === 'success' ? '\x1b[32m' : '\x1b[36m';
    const reset = '\x1b[0m';
    const metaStr = meta && Object.keys(meta).length > 0 ? ` ${JSON.stringify(meta)}` : '';
    return `${color}[${timestamp}] [${levelLabel}]${reset} ${message}${metaStr}`;
  }

  formatFileLine(level, message, meta) {
    const payload = {
      ts: new Date().toISOString(),
      level,
      msg: typeof message === 'string' ? message : String(message),
      ...(meta && Object.keys(meta).length > 0 ? { meta } : {})
    };
    return JSON.stringify(payload) + '\n';
  }

  writeToFile(filename, content) {
    if (!filename || !this.logsDir) return;
    
    try {
      fs.appendFileSync(filename, content);
    } catch (error) {
      console.error('Failed to write to log file:', error.message);
    }
  }

  log(level, message, meta = {}) {
    const normalizedLevel = level === 'success' ? 'success' : level;
    const safeMeta = this.safeSerializeMeta(meta) || {};
    const logFile = this.getLogFile();

    // File: JSONL for easy parsing/ingestion (respect LOG_LEVEL)
    if (this.shouldLog(normalizedLevel)) {
      this.writeToFile(logFile, this.formatFileLine(normalizedLevel, message, safeMeta));
    }

    // Console: human-friendly with colors; always log to console for development visibility
    const consoleMethod = normalizedLevel === 'error' ? 'error' : normalizedLevel === 'warn' ? 'warn' : 'log';
    console[consoleMethod](this.formatConsole(normalizedLevel, typeof message === 'string' ? message : String(message), safeMeta));
  }

  info(message, meta = {}) {
    this.log('info', message, meta);
  }

  success(message, meta = {}) {
    this.log('success', message, meta);
  }

  warn(message, meta = {}) {
    this.log('warn', message, meta);
  }

  error(message, meta = {}) {
    this.log('error', message, meta);
  }

  // Specific log types
  request(req, res, responseTime) {
    const meta = {
      method: req.method,
      url: req.url,
      status: res.statusCode,
      responseTime: `${responseTime}ms`,
      ip: req.ip || req.connection.remoteAddress,
      userAgent: req.get('User-Agent')
    };
    this.info(`${req.method} ${req.url}`, meta);
  }

  websocket(event, data = {}) {
    this.info(`WebSocket: ${event}`, data);
  }

  startup(message) {
    this.success(`🚀 ${message}`);
  }

  shutdown(message) {
    this.warn(`🛑 ${message}`);
  }
}

module.exports = FileLogger;
