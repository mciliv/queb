const fs = require('fs');
const path = require('path');

class FileLogger {
  constructor() {
    this.logsDir = path.join(__dirname, '..', '..', 'logs');
    this.ensureLogsDirectory();
  }

  ensureLogsDirectory() {
    if (!fs.existsSync(this.logsDir)) {
      fs.mkdirSync(this.logsDir, { recursive: true });
    }
  }

  getLogFile(type = 'server') {
    const today = new Date().toISOString().split('T')[0];
    return path.join(this.logsDir, `${type}-${today}.log`);
  }

  formatMessage(level, message, meta = {}) {
    const timestamp = new Date().toISOString();
    const metaStr = Object.keys(meta).length > 0 ? ` ${JSON.stringify(meta)}` : '';
    return `[${timestamp}] [${level.toUpperCase()}] ${message}${metaStr}\n`;
  }

  writeToFile(filename, content) {
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
