/**
 * Output Capture System for AI Command Integration
 * 
 * This module captures terminal output, logs, and error information
 * and provides it to AI commands for analysis and resolution.
 * It's completely separate from error handling and can be used
 * for any terminal output that needs AI analysis.
 * 
 * Key features:
 * - Capture any terminal output (stdout, stderr, logs)
 * - Format output for AI consumption
 * - Pass to any AI command or script
 * - Support multiple output sources
 * - Configurable output formatting
 */

const { spawn, exec } = require('child_process');
const { promisify } = require('util');
const fs = require('fs').promises;
const path = require('path');
const os = require('os');

const execAsync = promisify(exec);

class OutputCapture {
  constructor(config = {}) {
    this.logger = config.logger || console;
    this.tempDir = os.tmpdir();
    this.maxOutputSize = config.maxOutputSize || 1024 * 1024; // 1MB
    this.aiCommand = config.aiCommand || '~/Code/ai';
    this.autoTrigger = config.autoTrigger || false;
    this.triggerPatterns = config.triggerPatterns || [
      'error', 'failed', 'exception', 'timeout', 'crash',
      'EADDRINUSE', 'ENOTFOUND', 'ECONNREFUSED', 'EACCES',
      'TypeError', 'ReferenceError', 'SyntaxError'
    ];
    
    // Output buffers
    this.outputBuffer = [];
    this.errorBuffer = [];
    this.logBuffer = [];
    
    // Capture state
    this.isCapturing = false;
    this.captureStartTime = null;
  }

  /**
   * Start capturing output from a process
   */
  async captureProcess(command, args = [], options = {}) {
    this.logger.info(`🎬 Starting output capture for: ${command} ${args.join(' ')}`);
    
    this.isCapturing = true;
    this.captureStartTime = Date.now();
    this.outputBuffer = [];
    this.errorBuffer = [];
    this.logBuffer = [];
    
    return new Promise((resolve, reject) => {
      const process = spawn(command, args, {
        cwd: options.cwd || process.cwd(),
        env: { ...process.env, ...options.env },
        stdio: ['pipe', 'pipe', 'pipe']
      });
      
      let exitCode = null;
      let hasExited = false;
      
      // Capture stdout
      process.stdout.on('data', (data) => {
        const output = data.toString();
        this.outputBuffer.push({
          type: 'stdout',
          content: output,
          timestamp: new Date().toISOString()
        });
        this.logger.debug('STDOUT:', output.trim());
      });
      
      // Capture stderr
      process.stderr.on('data', (data) => {
        const output = data.toString();
        this.errorBuffer.push({
          type: 'stderr',
          content: output,
          timestamp: new Date().toISOString()
        });
        this.logger.debug('STDERR:', output.trim());
      });
      
      // Handle process exit
      process.on('exit', (code) => {
        exitCode = code;
        hasExited = true;
        this.isCapturing = false;
        
        const duration = Date.now() - this.captureStartTime;
        this.logger.info(`✅ Process completed with exit code ${code} (${duration}ms)`);
        
        resolve({
          exitCode,
          duration,
          output: this.outputBuffer,
          errors: this.errorBuffer,
          logs: this.logBuffer,
          success: code === 0
        });
      });
      
      // Handle process errors
      process.on('error', (error) => {
        if (!hasExited) {
          this.isCapturing = false;
          this.logger.error('❌ Process error:', error.message);
          reject(error);
        }
      });
      
      // Handle timeout
      if (options.timeout) {
        setTimeout(() => {
          if (!hasExited) {
            process.kill('SIGTERM');
            this.isCapturing = false;
            reject(new Error(`Process timeout after ${options.timeout}ms`));
          }
        }, options.timeout);
      }
    });
  }

  /**
   * Capture output from a command string
   */
  async captureCommand(command, options = {}) {
    this.logger.info(`🎬 Capturing command output: ${command}`);
    
    this.isCapturing = true;
    this.captureStartTime = Date.now();
    this.outputBuffer = [];
    this.errorBuffer = [];
    
    try {
      const { stdout, stderr } = await execAsync(command, {
        cwd: options.cwd || process.cwd(),
        timeout: options.timeout || 30000,
        maxBuffer: this.maxOutputSize
      });
      
      this.isCapturing = false;
      const duration = Date.now() - this.captureStartTime;
      
      // Process the output
      if (stdout) {
        this.outputBuffer.push({
          type: 'stdout',
          content: stdout,
          timestamp: new Date().toISOString()
        });
      }
      
      if (stderr) {
        this.errorBuffer.push({
          type: 'stderr',
          content: stderr,
          timestamp: new Date().toISOString()
        });
      }
      
      this.logger.info(`✅ Command completed (${duration}ms)`);
      
      return {
        exitCode: 0,
        duration,
        output: this.outputBuffer,
        errors: this.errorBuffer,
        logs: this.logBuffer,
        success: true
      };
      
    } catch (error) {
      this.isCapturing = false;
      const duration = Date.now() - this.captureStartTime;
      
      this.logger.error('❌ Command failed:', error.message);
      
      return {
        exitCode: error.code || 1,
        duration,
        output: this.outputBuffer,
        errors: this.errorBuffer,
        logs: this.logBuffer,
        success: false,
        error: error.message
      };
    }
  }

  /**
   * Add log entry to capture
   */
  addLog(level, message, context = {}) {
    this.logBuffer.push({
      type: 'log',
      level,
      content: message,
      context,
      timestamp: new Date().toISOString()
    });
    
    this.logger.debug(`[${level.toUpperCase()}] ${message}`);
  }

  /**
   * Check if captured output should trigger AI analysis
   */
  shouldTriggerAI(captureResult) {
    if (!this.autoTrigger) return false;
    
    const allOutput = [
      ...captureResult.output.map(o => o.content),
      ...captureResult.errors.map(e => e.content),
      ...captureResult.logs.map(l => l.content)
    ].join('\n').toLowerCase();
    
    return this.triggerPatterns.some(pattern => 
      allOutput.includes(pattern.toLowerCase())
    );
  }

  /**
   * Format captured output for AI consumption
   */
  formatForAI(captureResult, context = {}) {
    const timestamp = new Date().toISOString();
    const duration = captureResult.duration;
    const success = captureResult.success;
    
    let formatted = `# Terminal Output Analysis\n\n`;
    formatted += `**Timestamp:** ${timestamp}\n`;
    formatted += `**Duration:** ${duration}ms\n`;
    formatted += `**Success:** ${success ? 'Yes' : 'No'}\n`;
    formatted += `**Exit Code:** ${captureResult.exitCode}\n\n`;
    
    if (context.command) {
      formatted += `**Command:** ${context.command}\n`;
    }
    if (context.workingDirectory) {
      formatted += `**Working Directory:** ${context.workingDirectory}\n`;
    }
    if (context.environment) {
      formatted += `**Environment:** ${context.environment}\n`;
    }
    
    formatted += `\n## Output\n\n`;
    
    if (captureResult.output.length > 0) {
      formatted += `### Standard Output\n\`\`\`\n`;
      captureResult.output.forEach(entry => {
        formatted += `[${entry.timestamp}] ${entry.content}`;
      });
      formatted += `\n\`\`\`\n\n`;
    }
    
    if (captureResult.errors.length > 0) {
      formatted += `### Error Output\n\`\`\`\n`;
      captureResult.errors.forEach(entry => {
        formatted += `[${entry.timestamp}] ${entry.content}`;
      });
      formatted += `\n\`\`\`\n\n`;
    }
    
    if (captureResult.logs.length > 0) {
      formatted += `### Logs\n\`\`\`\n`;
      captureResult.logs.forEach(entry => {
        formatted += `[${entry.timestamp}] [${entry.level.toUpperCase()}] ${entry.content}`;
        if (entry.context && Object.keys(entry.context).length > 0) {
          formatted += `\nContext: ${JSON.stringify(entry.context, null, 2)}`;
        }
      });
      formatted += `\n\`\`\`\n\n`;
    }
    
    if (captureResult.error) {
      formatted += `### Process Error\n\`\`\`\n${captureResult.error}\n\`\`\`\n\n`;
    }
    
    formatted += `## Analysis Request\n\n`;
    formatted += `Please analyze this terminal output and provide:\n`;
    formatted += `1. What went wrong (if anything)\n`;
    formatted += `2. Root cause analysis\n`;
    formatted += `3. Specific fix recommendations\n`;
    formatted += `4. Commands to run to fix the issue\n`;
    formatted += `5. Prevention strategies\n`;
    
    return formatted;
  }

  /**
   * Pass captured output to AI command
   */
  async passToAI(captureResult, context = {}) {
    const formattedOutput = this.formatForAI(captureResult, context);
    
    // Write to temporary file
    const tempFile = path.join(this.tempDir, `output-capture-${Date.now()}.md`);
    await fs.writeFile(tempFile, formattedOutput);
    
    this.logger.info(`🤖 Passing output to AI: ${this.aiCommand}`);
    this.logger.debug('Formatted output:', formattedOutput);
    
    try {
      // Execute AI command with the formatted output
      const aiCommand = this.aiCommand.replace('~/Code/ai', path.join(os.homedir(), 'Code/ai'));
      const command = `cat "${tempFile}" | ${aiCommand}`;
      
      const { stdout, stderr } = await execAsync(command, {
        cwd: process.cwd(),
        timeout: 60000, // 1 minute timeout for AI
        maxBuffer: 1024 * 1024 * 10 // 10MB buffer
      });
      
      this.logger.info('✅ AI analysis completed');
      
      return {
        success: true,
        output: stdout,
        error: stderr,
        inputFile: tempFile,
        timestamp: new Date().toISOString()
      };
      
    } catch (error) {
      this.logger.error('❌ AI command failed:', error.message);
      
      return {
        success: false,
        error: error.message,
        inputFile: tempFile,
        timestamp: new Date().toISOString()
      };
    } finally {
      // Cleanup temp file
      try {
        await fs.unlink(tempFile);
      } catch (err) {
        this.logger.warn('Failed to cleanup temp file:', err.message);
      }
    }
  }

  /**
   * Capture and analyze in one step
   */
  async captureAndAnalyze(command, args = [], options = {}) {
    this.logger.info(`🎬 Capturing and analyzing: ${command} ${args.join(' ')}`);
    
    // Capture the output
    const captureResult = await this.captureProcess(command, args, options);
    
    // Check if we should trigger AI analysis
    if (this.shouldTriggerAI(captureResult) || options.forceAI) {
      this.logger.info('🤖 Triggering AI analysis...');
      
      const aiResult = await this.passToAI(captureResult, {
        command: `${command} ${args.join(' ')}`,
        workingDirectory: options.cwd || process.cwd(),
        environment: options.env ? Object.keys(options.env).join(', ') : 'default'
      });
      
      return {
        capture: captureResult,
        ai: aiResult
      };
    }
    
    return {
      capture: captureResult,
      ai: null
    };
  }

  /**
   * Get current capture state
   */
  getCaptureState() {
    return {
      isCapturing: this.isCapturing,
      startTime: this.captureStartTime,
      outputCount: this.outputBuffer.length,
      errorCount: this.errorBuffer.length,
      logCount: this.logBuffer.length
    };
  }

  /**
   * Clear all buffers
   */
  clearBuffers() {
    this.outputBuffer = [];
    this.errorBuffer = [];
    this.logBuffer = [];
    this.logger.debug('Cleared all capture buffers');
  }

  /**
   * Update configuration
   */
  updateConfig(newConfig) {
    Object.assign(this, newConfig);
    this.logger.info('Output capture configuration updated');
  }
}

module.exports = OutputCapture;