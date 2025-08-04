#!/usr/bin/env node

/**
 * File Watcher for Automated Testing
 * 
 * Watches files and automatically runs relevant tests when changes are detected
 * Usage: node test/watch/file-watcher.js [--patterns file1,file2]
 */

const chokidar = require('chokidar');
const path = require('path');
const { getTestsForFile } = require('../config/test-mapping');
const TestRunner = require('../config/test-runner');

class FileWatcher {
  constructor(options = {}) {
    this.testRunner = new TestRunner();
    this.isRunning = false;
    this.queuedFiles = new Set();
    this.debounceTimeout = null;
    this.debounceDelay = options.debounceDelay || 1000;
    this.patterns = options.patterns || [
      'backend/**/*.js',
      'frontend/**/*.js', 
      'frontend/**/*.html',
      'frontend/**/*.css',
      'chemistry/**/*.py',
      'package.json'
    ];
    
    console.log('🔍 File Watcher initialized');
    console.log(`📁 Watching patterns: ${this.patterns.join(', ')}`);
  }

  /**
   * Start watching files
   */
  start() {
    const watcher = chokidar.watch(this.patterns, {
      ignored: [
        '**/node_modules/**',
        '**/dist/**',
        '**/.git/**',
        '**/data/**',
        '**/logs/**',
        '**/*.log',
        '**/coverage/**'
      ],
      persistent: true,
      ignoreInitial: true
    });

    watcher.on('change', (filePath) => {
      this.handleFileChange(filePath);
    });

    watcher.on('add', (filePath) => {
      this.handleFileChange(filePath);
    });

    watcher.on('error', (error) => {
      console.error('❌ Watcher error:', error);
    });

    console.log('👀 Watching for file changes...');
    console.log('Press Ctrl+C to stop');

    // Handle graceful shutdown
    process.on('SIGINT', () => {
      console.log('\n🛑 Stopping file watcher...');
      watcher.close();
      this.testRunner.stopServer();
      process.exit(0);
    });

    return watcher;
  }

  /**
   * Handle file change event
   */
  handleFileChange(filePath) {
    const relativePath = path.relative(process.cwd(), filePath);
    console.log(`\n📝 File changed: ${relativePath}`);

    // Add to queue
    this.queuedFiles.add(relativePath);

    // Debounce to handle multiple rapid changes
    if (this.debounceTimeout) {
      clearTimeout(this.debounceTimeout);
    }

    this.debounceTimeout = setTimeout(() => {
      this.processQueuedChanges();
    }, this.debounceDelay);
  }

  /**
   * Process all queued file changes
   */
  async processQueuedChanges() {
    if (this.isRunning) {
      console.log('⏳ Tests already running, skipping...');
      return;
    }

    const changedFiles = Array.from(this.queuedFiles);
    this.queuedFiles.clear();

    if (changedFiles.length === 0) {
      return;
    }

    this.isRunning = true;

    try {
      console.log(`\n🚀 Processing ${changedFiles.length} changed file(s):`);
      changedFiles.forEach(file => console.log(`   • ${file}`));

      // Collect all tests to run
      let allTests = new Set();
      
      changedFiles.forEach(file => {
        const tests = getTestsForFile(file);
        console.log(`\n📁 ${file}:`);
        
        if (tests.length === 0) {
          console.log('   ⚠️ No tests mapped to this file');
        } else {
          console.log(`   🧪 ${tests.length} test suite(s):`);
          tests.forEach(test => {
            allTests.add(test);
            console.log(`      • ${test}`);
          });
        }
      });

      if (allTests.size === 0) {
        console.log('\n⏭️ No tests to run');
        return;
      }

      // Run tests
      const testList = Array.from(allTests);
      const success = await this.testRunner.runTests(testList);

      if (success) {
        console.log('\n✅ All tests passed! Watching for more changes...');
      } else {
        console.log('\n❌ Some tests failed. Fix issues and save files to re-run.');
      }

    } catch (error) {
      console.error('\n❌ Error running tests:', error.message);
    } finally {
      this.isRunning = false;
    }
  }

  /**
   * Run tests for specific files immediately
   */
  async runTestsForFiles(files) {
    this.isRunning = true;
    try {
      const success = await this.testRunner.runForFiles(files);
      return success;
    } finally {
      this.isRunning = false;
    }
  }
}

// CLI handling
async function main() {
  const args = process.argv.slice(2);

  if (args.includes('--help')) {
    console.log(`
File Watcher for Automated Testing

Usage:
  node test/watch/file-watcher.js [options]

Options:
  --patterns <file1,file2>  Custom file patterns to watch
  --debounce <ms>          Debounce delay in milliseconds (default: 1000)
  --test-files <file1,file2>  Run tests for specific files and exit
  --help                   Show this help

Examples:
  node test/watch/file-watcher.js
  node test/watch/file-watcher.js --patterns "backend/**/*.js,frontend/**/*.js"
  node test/watch/file-watcher.js --test-files backend/services/AtomPredictor.js
    `);
    return;
  }

  const patternsIndex = args.indexOf('--patterns');
  let patterns;
  if (patternsIndex !== -1 && args[patternsIndex + 1]) {
    patterns = args[patternsIndex + 1].split(',');
  }

  const debounceIndex = args.indexOf('--debounce');
  let debounceDelay;
  if (debounceIndex !== -1 && args[debounceIndex + 1]) {
    debounceDelay = parseInt(args[debounceIndex + 1], 10);
  }

  const testFilesIndex = args.indexOf('--test-files');
  if (testFilesIndex !== -1 && args[testFilesIndex + 1]) {
    const files = args[testFilesIndex + 1].split(',');
    const watcher = new FileWatcher({ patterns, debounceDelay });
    const success = await watcher.runTestsForFiles(files);
    process.exit(success ? 0 : 1);
    return;
  }

  // Start watching
  const watcher = new FileWatcher({ patterns, debounceDelay });
  watcher.start();
}

if (require.main === module) {
  main().catch(console.error);
}

module.exports = FileWatcher;