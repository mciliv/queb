const LogErrorPort = require('../ports/log-error-port');

/**
 * FILE LOGGER ADAPTER - INFRASTRUCTURE LAYER
 *
 * Adapts the FileLogger (infrastructure) to the LogErrorPort (domain contract).
 * Handles infrastructure concerns: file I/O, environment config, error handling.
 *
 * ISOLATION: Domain doesn't know about files, environment, or FileLogger specifics
 */
class FileLoggerAdapter extends LogErrorPort {
  constructor(fileLogger) {
    super();
    // INFRASTRUCTURE DEPENDENCY: Injected FileLogger instance
    this.fileLogger = fileLogger;
  }

  async log({ level, message, metadata }) {
    // INFRASTRUCTURE CONCERNS: File I/O, error handling, environment checks
    try {
      // ADAPT: Map domain log levels to FileLogger methods
      switch (level) {
        case 'error':
          this.fileLogger.error(message, metadata);
          break;
        case 'warn':
          this.fileLogger.warn(message, metadata);
          break;
        case 'info':
        default:
          this.fileLogger.info(message, metadata);
          break;
      }
    } catch (error) {
      // INFRASTRUCTURE ERROR HANDLING: Don't let logging failures crash the app
      console.error('Failed to write to log:', error.message);
      // Could fall back to console.log here if needed
    }
  }
}

module.exports = FileLoggerAdapter;