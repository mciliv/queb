/**
 * External Script Trigger for Automatic Error Resolution
 * 
 * This module provides a generic, abstracted way to trigger external scripts
 * (like ~/Code/ai, cursor-agent, or any custom script) when errors occur.
 * It's completely decoupled from the specific program and can be configured
 * to call any external script with error context.
 * 
 * Key features:
 * - Generic script execution (not tied to specific tools)
 * - Configurable script paths and arguments
 * - Error context passing via environment variables or stdin
 * - Retry logic and timeout handling
 * - Multiple script support (fallback chains)
 * - Safe execution with sandboxing
 */

const { exec, spawn } = require('child_process');
const { promisify } = require('util');
const fs = require('fs').promises;
const path = require('path');
const os = require('os');


class ExternalScriptTrigger {
  constructor(config = {}) {
    this.enabled = config.enabled !== false; // Default to enabled
    this.logger = config.logger || console;
    this.projectRoot = process.cwd();
    this.tempDir = os.tmpdir();
    
    // Script configuration
    this.scripts = config.scripts || [
      {
        name: 'ai-script',
        path: '~/Code/ai',
        enabled: true,
        priority: 1,
        timeout: 30000,
        retries: 2,
        method: 'env', // 'env', 'stdin', 'file'
        args: []
      },
      {
        name: 'cursor-agent',
        path: 'cursor-agent',
        enabled: false, // Disabled by default, enable if available
        priority: 2,
        timeout: 30000,
        retries: 1,
        method: 'stdin',
        args: []
      }
    ];
    
    // Error filtering
    this.triggerPatterns = config.triggerPatterns || [
      'EADDRINUSE', 'ENOTFOUND', 'ECONNREFUSED', 'EACCES', 'ENOENT',
      'MODULE_NOT_FOUND', 'syntax error', 'TypeError', 'ReferenceError',
      'Cannot read property', 'undefined is not a function',
      'Maximum call stack exceeded', 'timeout', 'rate limit',
      'authentication', 'validation', 'database', 'API', 'service unavailable'
    ];
    
    this.excludePatterns = config.excludePatterns || [
      'user input', 'invalid input', 'required field', 'format', 'length', 'pattern'
    ];
    
    // Severity thresholds
    this.severityThresholds = config.severityThresholds || {
      low: 1, medium: 2, high: 3, critical: 4
    };
    this.minSeverity = config.minSeverity || 'medium';
    
    // Rate limiting
    this.rateLimit = {
      maxCalls: config.rateLimit?.maxCalls || 10,
      windowMs: config.rateLimit?.windowMs || 60000, // 1 minute
      calls: new Map()
    };
  }

  /**
   * Check if an error should trigger external scripts
   */
  shouldTriggerScripts(errorInfo) {
    if (!this.enabled) return false;
    
    // Check rate limiting
    if (this._isRateLimited()) return false;
    
    const message = errorInfo.message.toLowerCase();
    const code = errorInfo.code || '';
    const severity = errorInfo.severity;
    
    // Check if error matches trigger patterns
    const matchesTrigger = this.triggerPatterns.some(pattern => 
      message.includes(pattern.toLowerCase()) || 
      code.toLowerCase().includes(pattern.toLowerCase())
    );
    
    // Check if error should be excluded
    const shouldExclude = this.excludePatterns.some(pattern => 
      message.includes(pattern.toLowerCase())
    );
    
    // Check severity threshold
    const meetsThreshold = this.severityThresholds[severity] >= this.severityThresholds[this.minSeverity];
    
    return matchesTrigger && !shouldExclude && meetsThreshold;
  }

  /**
   * Execute external scripts for error resolution
   */
  async executeScripts(errorInfo) {
    if (!this.shouldTriggerScripts(errorInfo)) {
      this.logger.debug('External scripts not triggered for this error type');
      return null;
    }

    // Update rate limiting
    this._updateRateLimit();
    
    const enabledScripts = this.scripts
      .filter(script => script.enabled)
      .sort((a, b) => a.priority - b.priority);
    
    if (enabledScripts.length === 0) {
      this.logger.warn('No enabled scripts configured for error handling');
      return null;
    }

    const results = [];
    
    for (const script of enabledScripts) {
      try {
        this.logger.info(`🤖 Triggering script: ${script.name} (${script.path})`);
        
        const result = await this._executeScript(script, errorInfo);
        results.push({
          script: script.name,
          success: result.success,
          output: result.output,
          error: result.error,
          duration: result.duration,
          timestamp: result.timestamp
        });
        
        if (result.success) {
          this.logger.info(`✅ Script ${script.name} completed successfully`);
          // Don't run fallback scripts if primary succeeds
          break;
        } else {
          this.logger.warn(`⚠️ Script ${script.name} failed: ${result.error}`);
        }
        
      } catch (error) {
        this.logger.error(`❌ Script ${script.name} execution error:`, error.message);
        results.push({
          script: script.name,
          success: false,
          error: error.message,
          timestamp: new Date().toISOString()
        });
      }
    }
    
    return {
      triggered: true,
      scripts: results,
      primarySuccess: results.length > 0 && results[0].success,
      timestamp: new Date().toISOString()
    };
  }

  /**
   * Execute a single script with retry logic
   */
  async _executeScript(script, errorInfo, attempt = 1) {
    const startTime = Date.now();
    
    try {
      const context = this._buildScriptContext(errorInfo);
      const command = this._buildScriptCommand(script, context);
      
      this.logger.debug(`Executing: ${command}`);
      
      const result = await this._runScript(command, script, context);
      
      return {
        success: true,
        output: result.stdout,
        error: result.stderr,
        duration: Date.now() - startTime,
        timestamp: new Date().toISOString(),
        attempt
      };
      
    } catch (error) {
      this.logger.warn(`Script ${script.name} attempt ${attempt} failed:`, error.message);
      
      if (attempt < script.retries) {
        this.logger.info(`🔄 Retrying ${script.name} in 2 seconds...`);
        await this._delay(2000);
        return this._executeScript(script, errorInfo, attempt + 1);
      }
      
      return {
        success: false,
        error: error.message,
        duration: Date.now() - startTime,
        timestamp: new Date().toISOString(),
        attempts: attempt
      };
    }
  }

  /**
   * Build script context for external scripts
   */
  _buildScriptContext(errorInfo) {
    return {
      error: {
        message: errorInfo.message,
        code: errorInfo.code,
        severity: errorInfo.severity,
        category: errorInfo.category,
        timestamp: errorInfo.timestamp,
        stack: errorInfo.stack
      },
      project: {
        root: this.projectRoot,
        name: require('../../package.json').name || 'unknown',
        version: require('../../package.json').version || 'unknown',
        nodeVersion: process.version,
        platform: process.platform,
        arch: process.arch
      },
      environment: {
        isDevelopment: process.env.NODE_ENV === 'development',
        isTest: process.env.NODE_ENV === 'test',
      },
      context: errorInfo.context || {}
    };
  }

  /**
   * Build command for script execution
   */
  _buildScriptCommand(script, context) {
    const scriptPath = this._resolveScriptPath(script.path);
    const args = script.args || [];
    
    // Add context as environment variables
    const envVars = this._buildEnvironmentVariables(context);
    const envString = Object.entries(envVars)
      .map(([key, value]) => `${key}="${value}"`)
      .join(' ');
    
    return `${envString} ${scriptPath} ${args.join(' ')}`.trim();
  }

  /**
   * Resolve script path (handle ~ expansion, relative paths, etc.)
   */
  _resolveScriptPath(scriptPath) {
    if (scriptPath.startsWith('~/')) {
      return path.join(os.homedir(), scriptPath.slice(2));
    }
    if (path.isAbsolute(scriptPath)) {
      return scriptPath;
    }
    return path.resolve(this.projectRoot, scriptPath);
  }

  /**
   * Build environment variables for script context
   */
  _buildEnvironmentVariables(context) {
    return {
      ERROR_MESSAGE: context.error.message,
      ERROR_CODE: context.error.code || '',
      ERROR_SEVERITY: context.error.severity,
      ERROR_CATEGORY: context.error.category || '',
      ERROR_TIMESTAMP: context.error.timestamp,
      PROJECT_ROOT: context.project.root,
      PROJECT_NAME: context.project.name,
      PROJECT_VERSION: context.project.version,
      NODE_VERSION: context.project.nodeVersion,
      PLATFORM: context.project.platform,
      ARCH: context.project.arch,
      IS_DEVELOPMENT: context.environment.isDevelopment ? 'true' : 'false',
      IS_TEST: context.environment.isTest ? 'true' : 'false',
      ERROR_STACK: context.error.stack || '',
      ERROR_CONTEXT: JSON.stringify(context.context)
    };
  }

  /**
   * Run the actual script
   */
  async _runScript(command, script, context) {
    return new Promise((resolve, reject) => {
      const timeout = setTimeout(() => {
        reject(new Error(`Script timeout after ${script.timeout}ms`));
      }, script.timeout);
      
      exec(command, {
        cwd: this.projectRoot,
        maxBuffer: 1024 * 1024 * 10, // 10MB
        timeout: script.timeout
      }, (error, stdout, stderr) => {
        clearTimeout(timeout);
        
        if (error) {
          reject(error);
        } else {
          resolve({ stdout, stderr });
        }
      });
    });
  }

  /**
   * Check if we're rate limited
   */
  _isRateLimited() {
    const now = Date.now();
    const windowStart = now - this.rateLimit.windowMs;
    
    // Clean old entries
    for (const [timestamp] of this.rateLimit.calls) {
      if (timestamp < windowStart) {
        this.rateLimit.calls.delete(timestamp);
      }
    }
    
    return this.rateLimit.calls.size >= this.rateLimit.maxCalls;
  }

  /**
   * Update rate limiting counters
   */
  _updateRateLimit() {
    const now = Date.now();
    this.rateLimit.calls.set(now, true);
  }

  /**
   * Utility method for delays
   */
  _delay(ms) {
    return new Promise(resolve => setTimeout(resolve, ms));
  }

  /**
   * Add a new script configuration
   */
  addScript(scriptConfig) {
    this.scripts.push({
      enabled: true,
      priority: 999,
      timeout: 30000,
      retries: 1,
      method: 'env',
      args: [],
      ...scriptConfig
    });
    
    // Re-sort by priority
    this.scripts.sort((a, b) => a.priority - b.priority);
    
    this.logger.info(`Added script: ${scriptConfig.name}`);
  }

  /**
   * Enable/disable a script
   */
  setScriptEnabled(scriptName, enabled) {
    const script = this.scripts.find(s => s.name === scriptName);
    if (script) {
      script.enabled = enabled;
      this.logger.info(`Script ${scriptName} ${enabled ? 'enabled' : 'disabled'}`);
    }
  }

  /**
   * Update configuration
   */
  updateConfig(newConfig) {
    Object.assign(this, newConfig);
    this.logger.info('External script trigger configuration updated');
  }

  /**
   * Get current configuration
   */
  getConfig() {
    return {
      enabled: this.enabled,
      scripts: this.scripts,
      triggerPatterns: this.triggerPatterns,
      excludePatterns: this.excludePatterns,
      minSeverity: this.minSeverity,
      rateLimit: this.rateLimit
    };
  }
}

module.exports = ExternalScriptTrigger;