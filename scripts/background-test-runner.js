#!/usr/bin/env node

const { spawn, exec } = require('child_process');
const fs = require('fs');
const path = require('path');

class BackgroundTestRunner {
  constructor() {
    this.isRunning = false;
    this.testQueue = new Set();
    this.lastResults = {};
  }

  log(message, type = 'info') {
    const timestamp = new Date().toLocaleTimeString();
    const symbols = {
      info: '🔍',
      success: '✅', 
      error: '❌',
      warning: '⚠️',
      test: '🧪'
    };
    console.log(`${symbols[type]} [${timestamp}] ${message}`);
  }

  async detectChanges() {
    return new Promise((resolve) => {
      exec('git diff --name-only HEAD~1 HEAD', (error, stdout) => {
        if (error) {
          // If no git history, check staged files
          exec('git diff --cached --name-only', (error2, stdout2) => {
            resolve(stdout2.split('\n').filter(f => f.trim()));
          });
        } else {
          resolve(stdout.split('\n').filter(f => f.trim()));
        }
      });
    });
  }

  analyzeBreakingChangeProbability(changedFiles) {
    const analysis = {
      probability: 'low',
      reasons: [],
      testSuites: new Set(['smoke'])
    };

    for (const file of changedFiles) {
      // High probability changes
      if (file.match(/backend\/(api|services|schemas)\/.*\.(js|ts)$/)) {
        analysis.probability = 'high';
        analysis.reasons.push(`API/Service change: ${file}`);
        analysis.testSuites.add('unit-backend');
        analysis.testSuites.add('integration');
      }
      
      if (file.match(/package\.json$|jest\.config|babel\.config/)) {
        analysis.probability = 'high';
        analysis.reasons.push(`Infrastructure change: ${file}`);
        analysis.testSuites.add('unit-backend');
        analysis.testSuites.add('unit-frontend');
      }

      // Medium probability changes
      if (file.match(/backend\/.*\.(js|ts)$/)) {
        analysis.probability = analysis.probability === 'high' ? 'high' : 'medium';
        analysis.reasons.push(`Backend change: ${file}`);
        analysis.testSuites.add('unit-backend');
      }

      if (file.match(/frontend\/.*\.(js|ts|html|css)$/)) {
        analysis.probability = analysis.probability === 'high' ? 'high' : 'medium';
        analysis.reasons.push(`Frontend change: ${file}`);
        analysis.testSuites.add('unit-frontend');
      }

      // Test file changes
      if (file.match(/test\/.*\.test\.(js|ts)$/)) {
        analysis.reasons.push(`Test change: ${file}`);
        if (file.includes('integration')) {
          analysis.testSuites.add('integration');
        }
      }
    }

    return analysis;
  }

  async runTestSuite(suiteName) {
    return new Promise((resolve) => {
      this.log(`Running ${suiteName} tests...`, 'test');
      
      const testCommands = {
        'smoke': 'npm run test:smoke',
        'unit-backend': 'npx jest --selectProjects unit-backend',
        'unit-frontend': 'npx jest --selectProjects unit-frontend', 
        'integration': 'npm run test:integration'
      };

      const cmd = testCommands[suiteName];
      if (!cmd) {
        this.log(`Unknown test suite: ${suiteName}`, 'error');
        resolve({ success: false, suite: suiteName });
        return;
      }

      const child = spawn('bash', ['-c', cmd], { 
        stdio: ['pipe', 'pipe', 'pipe'],
        env: { ...process.env, OPENAI_API_KEY: 'test-api-key-for-testing' }
      });

      let output = '';
      let errorOutput = '';

      child.stdout.on('data', (data) => {
        output += data.toString();
      });

      child.stderr.on('data', (data) => {
        errorOutput += data.toString();
      });

      child.on('close', (code) => {
        const success = code === 0;
        this.lastResults[suiteName] = {
          success,
          timestamp: new Date(),
          output,
          errorOutput
        };

        if (success) {
          this.log(`${suiteName} tests passed`, 'success');
        } else {
          this.log(`${suiteName} tests failed (exit code: ${code})`, 'error');
          if (errorOutput && !errorOutput.includes('mock issues')) {
            console.log('Error details:', errorOutput.slice(0, 500));
          }
        }

        resolve({ success, suite: suiteName, output, errorOutput });
      });
    });
  }

  async runIntelligentTests() {
    if (this.isRunning) {
      this.log('Tests already running, skipping...', 'warning');
      return;
    }

    this.isRunning = true;
    this.log('Starting intelligent background test runner...', 'info');

    try {
      const changedFiles = await this.detectChanges();
      
      if (changedFiles.length === 0) {
        this.log('No changes detected, running smoke tests only', 'info');
        await this.runTestSuite('smoke');
        return;
      }

      const analysis = this.analyzeBreakingChangeProbability(changedFiles);
      
      this.log(`Breaking change probability: ${analysis.probability.toUpperCase()}`, 'info');
      if (analysis.reasons.length > 0) {
        this.log(`Reasons: ${analysis.reasons.join(', ')}`, 'info');
      }

      const testSuites = Array.from(analysis.testSuites);
      this.log(`Running test suites: ${testSuites.join(', ')}`, 'test');

      // Run tests in order of importance
      const criticalSuites = ['smoke', 'unit-backend'];
      const optionalSuites = ['unit-frontend', 'integration'];

      // Run critical tests first
      for (const suite of criticalSuites) {
        if (testSuites.includes(suite)) {
          const result = await this.runTestSuite(suite);
          if (!result.success && suite === 'unit-backend') {
            this.log('Critical backend tests failed - stopping execution', 'error');
            break;
          }
        }
      }

      // Run optional tests (non-blocking)
      for (const suite of optionalSuites) {
        if (testSuites.includes(suite)) {
          await this.runTestSuite(suite);
        }
      }

      this.generateTestReport();

    } catch (error) {
      this.log(`Error running tests: ${error.message}`, 'error');
    } finally {
      this.isRunning = false;
    }
  }

  generateTestReport() {
    this.log('=== Test Results Summary ===', 'info');
    
    let allPassed = true;
    for (const [suite, result] of Object.entries(this.lastResults)) {
      const status = result.success ? '✅ PASS' : '❌ FAIL';
      this.log(`${suite}: ${status}`, result.success ? 'success' : 'error');
      if (!result.success && !suite.includes('frontend')) {
        allPassed = false;
      }
    }

    if (allPassed) {
      this.log('All critical tests passing! 🎉', 'success');
    } else {
      this.log('Some critical tests failed - check output above', 'warning');
    }
  }

  startWatcher() {
    this.log('Starting file watcher for continuous testing...', 'info');
    
    // Watch for git commits
    const gitWatcher = fs.watch('.git/refs/heads', { recursive: true }, (eventType, filename) => {
      if (eventType === 'change') {
        this.log('Git change detected, running tests...', 'info');
        setTimeout(() => this.runIntelligentTests(), 1000);
      }
    });

    // Watch for staged changes
    setInterval(() => {
      exec('git diff --cached --name-only', (error, stdout) => {
        if (!error && stdout.trim() && !this.isRunning) {
          const stagedFiles = stdout.split('\n').filter(f => f.trim());
          if (stagedFiles.length > 0) {
            this.log(`Staged files detected: ${stagedFiles.length} files`, 'info');
          }
        }
      });
    }, 30000); // Check every 30 seconds
  }
}

// CLI interface
if (require.main === module) {
  const runner = new BackgroundTestRunner();
  
  const command = process.argv[2];
  
  switch (command) {
    case 'run':
      runner.runIntelligentTests();
      break;
    case 'watch':
      runner.startWatcher();
      break;
    case 'status':
      runner.generateTestReport();
      break;
    default:
      console.log('Usage: node background-test-runner.js [run|watch|status]');
      console.log('  run    - Run tests based on changes');
      console.log('  watch  - Start continuous testing');
      console.log('  status - Show last test results');
  }
}

module.exports = BackgroundTestRunner;