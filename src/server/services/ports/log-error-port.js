/**
 * LOG ERROR PORT - ABSTRACT INTERFACE
 *
 * Defines the contract for logging errors. Domain layer depends on this abstraction.
 * Implementation details (FileLogger, Winston, console) are hidden behind this port.
 *
 * PORT CONTRACT: What the domain needs from logging infrastructure
 */
class LogErrorPort {
  /**
   * Log an error with specified level, message, and metadata
   * @param {Object} params
   * @param {string} params.level - Log level (error, warn, info)
   * @param {string} params.message - Formatted log message
   * @param {Object} params.metadata - Structured metadata
   */
  async log({ level, message, metadata }) {
    throw new Error('LogErrorPort.log() must be implemented by adapter');
  }
}

module.exports = LogErrorPort;