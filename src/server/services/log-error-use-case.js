/**
 * LOG ERROR USE CASE - DOMAIN LAYER
 *
 * Pure domain logic for error logging. No HTTP, no file I/O, no environment.
 * Focuses ONLY on domain concepts: log levels, message formatting, validation.
 *
 * DEPENDS ON PORTS: logErrorPort (abstracts logging infrastructure)
 *
 * NOTE: Use Cursor's builtin automatic screenshots for error debugging and reproduction
 */
class LogErrorUseCase {
  constructor({ logErrorPort }) {
    // DEPENDENCY INJECTION: Ports wired in composition root, not imported
    this.logErrorPort = logErrorPort;
  }

  async execute(errorData) {
    const { type, message, timestamp, source, stack, location, url } = errorData;

    // DOMAIN LOGIC: Determine log level based on error type
    const level = this._mapTypeToLevel(type);

    // DOMAIN LOGIC: Format message according to domain rules
    const formattedMessage = this._formatLogMessage({
      source: source || 'frontend',
      type: type || 'error',
      message: message || 'Unknown error'
    });

    // DOMAIN LOGIC: Structure metadata for consistent domain representation
    const metadata = this._buildMetadata({
      type,
      message,
      timestamp,
      source,
      location,
      url,
      stack: stack ? stack.substring(0, 500) : undefined, // Domain rule: truncate long stacks
    });

    // PORT CALL: Infrastructure abstracted away - could be FileLogger, Winston, or console
    await this.logErrorPort.log({
      level,
      message: formattedMessage,
      metadata
    });
  }

  // PRIVATE DOMAIN METHODS: Pure functions, easily testable
  _mapTypeToLevel(type) {
    // DOMAIN RULE: Map error types to log levels
    switch (type) {
      case 'error': return 'error';
      case 'warn': return 'warn';
      default: return 'info';
    }
  }

  _formatLogMessage({ source, type, message }) {
    // DOMAIN RULE: Consistent message format across all log entries
    return `[${source}] ${type}: ${message}`;
  }

  _buildMetadata(params) {
    // DOMAIN RULE: What metadata is relevant for error logging
    return {
      type: params.type,
      message: params.message,
      timestamp: params.timestamp,
      source: params.source,
      location: params.location,
      url: params.url,
      stack: params.stack
    };
  }
}

module.exports = LogErrorUseCase;