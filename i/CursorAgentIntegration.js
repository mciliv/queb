/**
 * Cursor Agent Integration for Automatic Error Resolution
 *
 * This module provides automatic cursor-agent integration for error scenarios.
 * It extends ExternalScriptTrigger to add cursor-agent-specific intelligent
 * prompt generation while leveraging the unified script execution infrastructure.
 *
 * Key features:
 * - Intelligent error prompt generation for cursor-agent
 * - Context-aware error resolution suggestions
 * - Unified script execution via ExternalScriptTrigger base class
 */

const ExternalScriptTrigger = require('./ExternalScriptTrigger');
const services = require('./services');

class CursorAgentIntegration extends ExternalScriptTrigger {
  constructor(config = {}) {
    // Configure cursor-agent as the primary script
    const cursorAgentConfig = {
      enabled: config.enabled !== false,
      logger: config.logger || console,
      minSeverity: config.autoFixThreshold || 'medium',
      scripts: [
        {
          name: 'cursor-agent',
          path: config.cursorAgentPath || 'cursor-agent',
          enabled: true,
          priority: 1,
          timeout: config.timeout || 30000,
          retries: config.maxRetries || 3,
          method: 'args', // Pass prompt as argument
          args: [] // Will be dynamically generated
        }
      ],
      ...config
    };

    super(cursorAgentConfig);

    this.projectRoot = process.cwd();
    this.retryDelay = config.retryDelay || 2000;
  }

  /**
   * Build intelligent prompt for cursor-agent
   */
  buildErrorPrompt(errorInfo, context) {
    const severity = errorInfo.severity.toUpperCase();
    const category = errorInfo.category || 'unknown';

    return `I'm experiencing a ${severity} severity ${category} error in my Node.js application.

ERROR DETAILS:
- Message: ${errorInfo.message}
- Code: ${errorInfo.code || 'N/A'}
- Severity: ${errorInfo.severity}
- Category: ${category}
- Timestamp: ${errorInfo.timestamp}

PROJECT CONTEXT:
- Project Root: ${context.project.root}
- Node Version: ${context.project.nodeVersion}
- Platform: ${context.project.platform}
- Environment: ${context.environment.isDevelopment ? 'Development' : 'Production'}

${errorInfo.stack ? `STACK TRACE:\n${errorInfo.stack}` : ''}

Please analyze this error and provide:
1. Root cause analysis
2. Specific fix recommendations
3. Code changes needed (if any)
4. Prevention strategies
5. Alternative solutions if the primary fix doesn't work

Focus on actionable solutions that can be implemented immediately.`;
  }

  /**
   * Override script command building to inject intelligent prompt
   */
  _buildScriptCommand(script, context) {
    if (script.name === 'cursor-agent') {
      // Build cursor-agent specific prompt
      const errorInfo = {
        message: context.error.message,
        code: context.error.code,
        severity: context.error.severity,
        category: context.error.category,
        timestamp: context.error.timestamp,
        stack: context.error.stack,
        context: context.context
      };

      const prompt = this.buildErrorPrompt(errorInfo, context);
      const scriptPath = this._resolveScriptPath(script.path);

      // Escape the prompt for shell execution
      const escapedPrompt = prompt.replace(/'/g, "'\"'\"'");

      return `${scriptPath} "${escapedPrompt}"`;
    }

    // Fall back to parent implementation for other scripts
    return super._buildScriptCommand(script, context);
  }

  /**
   * Execute cursor-agent (delegates to base class executeScripts)
   */
  async executeCursorAgent(errorInfo, attempt = 1) {
    return this.executeScripts(errorInfo);
  }

  /**
   * Check if cursor-agent should be triggered (delegates to base class)
   */
  shouldTriggerCursorAgent(errorInfo) {
    return this.shouldTriggerScripts(errorInfo);
  }

  /**
   * Process cursor-agent response
   */
  processCursorAgentResponse(response) {
    if (!response || !response.primarySuccess) {
      return {
        success: false,
        error: response?.scripts?.[0]?.error || 'No response from cursor-agent',
        timestamp: response?.timestamp || new Date().toISOString()
      };
    }

    return {
      success: true,
      output: response.scripts[0].output,
      timestamp: response.timestamp
    };
  }

  /**
   * Override delay to use custom retry delay
   */
  _delay(ms) {
    return super._delay(this.retryDelay);
  }
}

module.exports = CursorAgentIntegration;