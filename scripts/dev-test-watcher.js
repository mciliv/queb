#!/usr/bin/env node

const fs = require('fs');
const path = require('path');
const { exec } = require('child_process');
const BackgroundTestRunner = require('./background-test-runner');

class DevTestWatcher {
  constructor() {
    this.testRunner = new BackgroundTestRunner();
    this.watchedDirs = [
      'backend',
      'frontend', 
      'test',
      'package.json',
      'jest.config.js'
    ];
    this.debounceTimer = null;
    this.isRunning = false;
  }

  log(message, type = 'info') {
    const timestamp = new Date().toLocaleTimeString();
    const symbols = {
      info: '👀',
      success: '✅', 
      error: '❌',
      warning: '⚠️',
      change: '📝'
    };
    console.log(`${symbols[type]} [${timestamp}] ${message}`);
  }

  debounceTest(changedFile) {
    clearTimeout(this.debounceTimer);
    
    this.debounceTimer = setTimeout(() => {
      if (!this.isRunning) {
        this.runTestsForFile(changedFile);
      }
    }, 2000); // Wait 2 seconds for multiple changes
  }

  async runTestsForFile(changedFile) {
    this.isRunning = true;
    this.log(`Change detected: ${changedFile}`, 'change');

    // Quick smoke test for immediate feedback
    if (changedFile.includes('backend') || changedFile.includes('test')) {
      this.log('Running quick smoke test...', 'info');
      
      try {
        await new Promise((resolve) => {
          exec('npm run test:smoke', (error, stdout, stderr) => {
            if (error) {
              this.log('Smoke test failed!', 'error');
              console.log(stderr);
            } else {
              this.log('Smoke test passed ✅', 'success');
            }
            resolve();
          });
        });
        
        // Run background tests
        setTimeout(() => {
          this.testRunner.runIntelligentTests();
        }, 1000);
        
      } catch (error) {
        this.log(`Error running tests: ${error.message}`, 'error');
      }
    }
    
    this.isRunning = false;
  }

  startWatching() {
    this.log('🚀 Starting development test watcher...', 'info');
    this.log(`Watching directories: ${this.watchedDirs.join(', ')}`, 'info');
    this.log('Press Ctrl+C to stop', 'info');
    console.log();

    const watchers = [];

    for (const dir of this.watchedDirs) {
      if (fs.existsSync(dir)) {
        const stats = fs.statSync(dir);
        
        if (stats.isDirectory()) {
          const watcher = fs.watch(dir, { recursive: true }, (eventType, filename) => {
            if (filename && (filename.endsWith('.js') || filename.endsWith('.ts') || filename.endsWith('.json'))) {
              this.debounceTest(path.join(dir, filename));
            }
          });
          watchers.push(watcher);
        } else {
          // Single file
          const watcher = fs.watch(dir, (eventType, filename) => {
            this.debounceTest(dir);
          });
          watchers.push(watcher);
        }
      }
    }

    // Graceful shutdown
    process.on('SIGINT', () => {
      this.log('Shutting down test watcher...', 'info');
      watchers.forEach(watcher => watcher.close());
      process.exit(0);
    });

    // Status updates every 30 seconds
    setInterval(() => {
      if (!this.isRunning) {
        this.log('Watching for changes... (all quiet)', 'info');
      }
    }, 30000);
  }

  displayStatus() {
    this.log('=== Development Test Watcher Status ===', 'info');
    this.log(`Watched directories: ${this.watchedDirs.length}`, 'info');
    this.log(`Currently running: ${this.isRunning ? 'Yes' : 'No'}`, 'info');
    
    if (this.testRunner.lastResults && Object.keys(this.testRunner.lastResults).length > 0) {
      this.log('Last test results:', 'info');
      this.testRunner.generateTestReport();
    } else {
      this.log('No test results yet', 'info');
    }
  }
}

// CLI interface
if (require.main === module) {
  const watcher = new DevTestWatcher();
  
  const command = process.argv[2];
  
  switch (command) {
    case 'start':
      watcher.startWatching();
      break;
    case 'status':
      watcher.displayStatus();
      break;
    default:
      console.log('Usage: node dev-test-watcher.js [start|status]');
      console.log('  start  - Start watching files for changes');
      console.log('  status - Show watcher status');
  }
}

module.exports = DevTestWatcher;