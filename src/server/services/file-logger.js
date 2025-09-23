const fs = require('fs');
const path = require('path');

class FileLogger {
  constructor() {
    // Use a more robust logs directory path that works in serverless environments
    // Always use the project root for logs, not the current working directory
    const projectRoot = this.findProjectRoot();
    this.logsDir = process.env.LOGS_DIR || path.join(projectRoot, 'logs');
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
    
    // Fallback to process.cwd() if we can't find package.json
    return process.cwd();
  }

  ensureLogsDirectory() {
    try {
      if (!fs.existsSync(this.logsDir)) {
        fs.mkdirSync(this.logsDir, { recursive: true });
      }
    } catch (error) {
      // In serverless environments, fall back to console logging only
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
        this.writeToFile(this.getLogFile(), this.formatMessage('INFO', `🧹 Cleared ${clearedCount} old log file(s) from previous session`));
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

  formatMessage(level, message, meta = {}) {
    const timestamp = new Date().toISOString();
    const metaStr = Object.keys(meta).length > 0 ? ` ${JSON.stringify(meta)}` : '';
    return `[${timestamp}] [${level.toUpperCase()}] ${message}${metaStr}\n`;
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
    const formattedMessage = this.formatMessage(level, message, meta);
    const logFile = this.getLogFile();

    // Write to file
    this.writeToFile(logFile, formattedMessage);

    // Also write to console for development
    const consoleMethod = level === 'error' ? 'error' :
                         level === 'warn' ? 'warn' : 'log';
    console[consoleMethod](formattedMessage.trim());
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

module.exports = new FileLogger();
