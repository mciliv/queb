#!/usr/bin/env node
/**
 * Node.js wrapper for running commands with AI analysis
 * 
 * This script can be used to run any command and automatically
 * pass the output to an AI command for analysis when errors occur.
 * 
 * Usage:
 *   node scripts/run-with-ai.js "npm start"
 *   node scripts/run-with-ai.js "node index.js" --force-ai
 *   node scripts/run-with-ai.js "npm test" --ai-command "~/Code/ai"
 */

const { spawn, exec } = require('child_process');
const { promisify } = require('util');
const fs = require('fs').promises;
const path = require('path');
const os = require('os');

const execAsync = promisify(exec);

class CommandRunner {
  constructor(options = {}) {
    this.aiCommand = options.aiCommand || '~/Code/ai';
    this.autoTrigger = options.autoTrigger !== false;
    this.forceAI = options.forceAI || false;
    this.timeout = options.timeout || 300000; // 5 minutes
    this.logger = options.logger || console;
    
    this.triggerPatterns = [
      'error', 'failed', 'exception', 'timeout', 'crash',
      'EADDRINUSE', 'ENOTFOUND', 'ECONNREFUSED', 'EACCES',
      'TypeError', 'ReferenceError', 'SyntaxError',
      'Cannot read property', 'undefined is not a function',
      'Maximum call stack exceeded', 'RangeError',
      'listen EADDRINUSE', 'address already in use',
      'Module not found', 'Cannot find module'
    ];
  }

  /**
   * Run a command and capture its output
   */
  async runCommand(command, args = []) {
    this.logger.info(`🎬 Running: ${command} ${args.join(' ')}`);
    
    const startTime = Date.now();
    const output = [];
    const errors = [];
    
    return new Promise((resolve, reject) => {
      const process = spawn(command, args, {
        stdio: ['pipe', 'pipe', 'pipe'],
        cwd: process.cwd()
      });
      
      let hasExited = false;
      
      // Capture stdout
      process.stdout.on('data', (data) => {
        const line = data.toString();
        output.push({ type: 'stdout', content: line, timestamp: new Date() });
        process.stdout.write(line); // Also show in real-time
      });
      
      // Capture stderr
      process.stderr.on('data', (data) => {
        const line = data.toString();
        errors.push({ type: 'stderr', content: line, timestamp: new Date() });
        process.stderr.write(line); // Also show in real-time
      });
      
      // Handle process exit
      process.on('exit', (code) => {
        if (!hasExited) {
          hasExited = true;
          const duration = Date.now() - startTime;
          
          resolve({
            exitCode: code,
            duration,
            output,
            errors,
            success: code === 0,
            command: `${command} ${args.join(' ')}`
          });
        }
      });
      
      // Handle process errors
      process.on('error', (error) => {
        if (!hasExited) {
          hasExited = true;
          reject(error);
        }
      });
      
      // Handle timeout
      setTimeout(() => {
        if (!hasExited) {
          process.kill('SIGTERM');
          hasExited = true;
          reject(new Error(`Command timeout after ${this.timeout}ms`));
        }
      }, this.timeout);
    });
  }

  /**
   * Check if output should trigger AI analysis
   */
  shouldTriggerAI(result) {
    if (this.forceAI) return true;
    if (!this.autoTrigger) return false;
    if (result.success) return false;
    
    const allOutput = [
      ...result.output.map(o => o.content),
      ...result.errors.map(e => e.content)
    ].join('\n').toLowerCase();
    
    return this.triggerPatterns.some(pattern => 
      allOutput.includes(pattern.toLowerCase())
    );
  }

  /**
   * Format output for AI consumption
   */
  formatForAI(result) {
    const timestamp = new Date().toISOString();
    const duration = result.duration;
    const success = result.success;
    
    let formatted = `# Terminal Output Analysis\n\n`;
    formatted += `**Command:** ${result.command}\n`;
    formatted += `**Exit Code:** ${result.exitCode}\n`;
    formatted += `**Duration:** ${duration}ms\n`;
    formatted += `**Success:** ${success ? 'Yes' : 'No'}\n`;
    formatted += `**Timestamp:** ${timestamp}\n\n`;
    
    formatted += `## Standard Output\n\`\`\`\n`;
    result.output.forEach(entry => {
      formatted += `[${entry.timestamp.toISOString()}] ${entry.content}`;
    });
    formatted += `\n\`\`\`\n\n`;
    
    formatted += `## Error Output\n\`\`\`\n`;
    result.errors.forEach(entry => {
      formatted += `[${entry.timestamp.toISOString()}] ${entry.content}`;
    });
    formatted += `\n\`\`\`\n\n`;
    
    formatted += `## Analysis Request\n\n`;
    formatted += `Please analyze this terminal output and provide:\n`;
    formatted += `1. What went wrong (if anything)\n`;
    formatted += `2. Root cause analysis\n`;
    formatted += `3. Specific fix recommendations\n`;
    formatted += `4. Commands to run to fix the issue\n`;
    formatted += `5. Prevention strategies\n\n`;
    formatted += `Focus on actionable solutions that can be implemented immediately.`;
    
    return formatted;
  }

  /**
   * Pass output to AI command
   */
  async passToAI(result) {
    const formattedOutput = this.formatForAI(result);
    
    // Write to temporary file
    const tempFile = path.join(os.tmpdir(), `output-capture-${Date.now()}.md`);
    await fs.writeFile(tempFile, formattedOutput);
    
    this.logger.info(`🤖 Passing output to AI: ${this.aiCommand}`);
    
    try {
      // Resolve AI command path
      const aiCommand = this.aiCommand.replace('~/Code/ai', path.join(os.homedir(), 'Code/ai'));
      
      // Execute AI command
      const { stdout, stderr } = await execAsync(`cat "${tempFile}" | ${aiCommand}`, {
        cwd: process.cwd(),
        timeout: 60000, // 1 minute timeout for AI
        maxBuffer: 1024 * 1024 * 10 // 10MB buffer
      });
      
      this.logger.info('✅ AI analysis completed');
      
      return {
        success: true,
        output: stdout,
        error: stderr,
        inputFile: tempFile
      };
      
    } catch (error) {
      this.logger.error('❌ AI command failed:', error.message);
      
      return {
        success: false,
        error: error.message,
        inputFile: tempFile
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
   * Run command and analyze with AI if needed
   */
  async runWithAI(command, args = []) {
    try {
      // Run the command
      const result = await this.runCommand(command, args);
      
      // Check if we should trigger AI analysis
      if (this.shouldTriggerAI(result)) {
        this.logger.info('🤖 Triggering AI analysis...');
        
        const aiResult = await this.passToAI(result);
        
        return {
          command: result,
          ai: aiResult
        };
      } else {
        this.logger.info('✅ Command completed successfully');
        return {
          command: result,
          ai: null
        };
      }
      
    } catch (error) {
      this.logger.error('❌ Command execution failed:', error.message);
      throw error;
    }
  }
}

// CLI interface
async function main() {
  const args = process.argv.slice(2);
  
  if (args.length === 0) {
    console.log('Usage: node scripts/run-with-ai.js "command" [args...] [options]');
    console.log('');
    console.log('Options:');
    console.log('  --force-ai              Force AI analysis even on success');
    console.log('  --no-auto-trigger       Disable automatic AI triggering');
    console.log('  --ai-command <path>     Specify AI command path');
    console.log('  --timeout <ms>          Set command timeout in milliseconds');
    console.log('');
    console.log('Examples:');
    console.log('  node scripts/run-with-ai.js "npm start"');
    console.log('  node scripts/run-with-ai.js "node index.js" --force-ai');
    console.log('  node scripts/run-with-ai.js "npm test" --ai-command "~/Code/ai"');
    process.exit(1);
  }
  
  // Parse options
  const options = {
    forceAI: false,
    autoTrigger: true,
    aiCommand: '~/Code/ai',
    timeout: 300000
  };
  
  const commandArgs = [];
  
  for (let i = 0; i < args.length; i++) {
    const arg = args[i];
    
    switch (arg) {
      case '--force-ai':
        options.forceAI = true;
        break;
      case '--no-auto-trigger':
        options.autoTrigger = false;
        break;
      case '--ai-command':
        options.aiCommand = args[++i];
        break;
      case '--timeout':
        options.timeout = parseInt(args[++i]);
        break;
      default:
        commandArgs.push(arg);
        break;
    }
  }
  
  if (commandArgs.length === 0) {
    console.error('Error: No command specified');
    process.exit(1);
  }
  
  // Create runner and execute
  const runner = new CommandRunner(options);
  
  try {
    const [command, ...args] = commandArgs;
    const result = await runner.runWithAI(command, args);
    
    // Exit with the same code as the original command
    process.exit(result.command.exitCode);
    
  } catch (error) {
    console.error('Fatal error:', error.message);
    process.exit(1);
  }
}

// Run if called directly
if (require.main === module) {
  main();
}

module.exports = CommandRunner;