#!/usr/bin/env node

/**
 * Smart Test Runner
 * 
 * Runs tests based on file changes, respecting dependency layers
 * Usage:
 *   node test/config/test-runner.js --files backend/services/AtomPredictor.js
 *   node test/config/test-runner.js --layer unit
 *   node test/config/test-runner.js --all
 *   node test/config/test-runner.js --integration-setup
 */

const { execSync, spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const { 
  getTestsForFile, 
  getTestsForLayer, 
  TEST_LAYERS,
  LAYER_TEST_SUITES 
} = require('./test-mapping');

class TestRunner {
  constructor() {
    this.serverProcess = null;
    this.isVerbose = process.argv.includes('--verbose');
  }

  log(message) {
    if (this.isVerbose) {
      console.log(`[TestRunner] ${message}`);
    }
  }

  /**
   * Start development server for integration tests
   */
  async startServer() {
    return new Promise((resolve, reject) => {
      this.log('Starting development server for integration tests...');
      
      // Kill any existing processes on port 8080
      try {
        execSync('lsof -ti:8080 | xargs kill -9', { stdio: 'ignore' });
        this.log('Killed existing processes on port 8080');
      } catch (e) {
        // No processes to kill
      }

      // Start server in background
      this.serverProcess = spawn('node', ['backend/api/server.js'], {
        env: { ...process.env, NODE_ENV: 'test', PORT: '8080' },
        stdio: this.isVerbose ? 'inherit' : 'ignore'
      });

      // Wait for server to be ready
      const checkServer = () => {
        try {
          execSync('curl -f http://localhost:8080 > /dev/null 2>&1');
          this.log('Server is ready');
          resolve();
        } catch (e) {
          setTimeout(checkServer, 500);
        }
      };

      setTimeout(checkServer, 2000);

      // Handle process errors
      this.serverProcess.on('error', (err) => {
        this.log(`Server error: ${err.message}`);
        reject(err);
      });
    });
  }

  /**
   * Stop development server
   */
  stopServer() {
    if (this.serverProcess) {
      this.log('Stopping development server...');
      this.serverProcess.kill();
      this.serverProcess = null;
    }
  }

  /**
   * Run specific test files
   */
  async runTests(testFiles, options = {}) {
    if (testFiles.length === 0) {
      console.log('No tests to run');
      return true;
    }

    console.log(`\n🧪 Running ${testFiles.length} test suite(s):`);
    testFiles.forEach(file => console.log(`  • ${file}`));

    const needsServer = testFiles.some(file => 
      file.includes('integration') && !file.includes('molecular-accuracy')
    );

    try {
      if (needsServer && !options.skipServer) {
        await this.startServer();
      }

      // Filter out non-existent files
      const existingTests = testFiles.filter(file => {
        const fullPath = path.join(process.cwd(), file);
        return fs.existsSync(fullPath);
      });

      if (existingTests.length === 0) {
        console.log('⚠️ No test files found');
        return true;
      }

      const jestArgs = [
        '--config', 'test/jest.config.js',
        '--testPathPattern', `"(${existingTests.join('|')})"`,
        '--runInBand', // Run tests serially for integration tests
      ];

      if (this.isVerbose) {
        jestArgs.push('--verbose');
      }

      console.log(`\n▶️ npx jest ${jestArgs.join(' ')}\n`);

      const result = execSync(`npx jest ${jestArgs.join(' ')}`, {
        stdio: 'inherit',
        encoding: 'utf8'
      });

      console.log('\n✅ All tests passed!');
      return true;

    } catch (error) {
      console.log(`\n❌ Tests failed with exit code: ${error.status}`);
      return false;
    } finally {
      if (needsServer && !options.skipServer) {
        this.stopServer();
      }
    }
  }

  /**
   * Run tests for changed files
   */
  async runForFiles(files) {
    let allTests = new Set();

    files.forEach(file => {
      const tests = getTestsForFile(file);
      console.log(`\n📁 ${file} → ${tests.length} test suite(s)`);
      tests.forEach(test => {
        // Only add existing files
        if (fs.existsSync(test)) {
          allTests.add(test);
          console.log(`   ✓ ${test}`);
        } else {
          console.log(`   ⚠️ ${test} (not found)`);
        }
      });
    });

    if (allTests.size === 0) {
      console.log('\n⚠️ No existing test files found');
      return true;
    }

    return this.runTests(Array.from(allTests));
  }

  /**
   * Run tests for a specific layer
   */
  async runForLayer(layer) {
    const tests = getTestsForLayer(layer);
    console.log(`\n🏗️ Running ${layer} layer tests (${tests.length} suites)`);
    return this.runTests(tests);
  }

  /**
   * Run all tests in dependency order
   */
  async runAll() {
    console.log('\n🎯 Running all tests in dependency order...');
    
    const layers = Object.values(TEST_LAYERS);
    let allPassed = true;

    for (const layer of layers) {
      console.log(`\n🏗️ === ${layer.toUpperCase()} LAYER ===`);
      const layerTests = LAYER_TEST_SUITES[layer];
      
      if (layerTests.length === 0) {
        console.log(`⏭️ No tests in ${layer} layer`);
        continue;
      }

      const passed = await this.runTests(layerTests, { 
        skipServer: layer !== TEST_LAYERS.INTEGRATION 
      });
      
      if (!passed) {
        console.log(`❌ ${layer} layer tests failed - stopping execution`);
        allPassed = false;
        break;
      }
    }

    return allPassed;
  }

  /**
   * Setup integration test environment
   */
  async setupIntegration() {
    console.log('🔧 Setting up integration test environment...');
    
    try {
      await this.startServer();
      console.log('✅ Integration environment ready');
      console.log('🌐 Server running at http://localhost:8080');
      console.log('Press Ctrl+C to stop');
      
      // Keep running until interrupted
      process.on('SIGINT', () => {
        console.log('\n🛑 Shutting down...');
        this.stopServer();
        process.exit(0);
      });

    } catch (error) {
      console.error('❌ Failed to setup integration environment:', error.message);
      this.stopServer();
      process.exit(1);
    }
  }
}

// CLI handling
async function main() {
  const runner = new TestRunner();
  const args = process.argv.slice(2);

  try {
    if (args.includes('--help')) {
      console.log(`
Smart Test Runner

Usage:
  node test/config/test-runner.js [options]

Options:
  --files <file1,file2>     Run tests for specific files
  --layer <layer>           Run tests for specific layer (${Object.values(TEST_LAYERS).join(', ')})
  --all                     Run all tests in dependency order
  --integration-setup       Start server for manual integration testing
  --verbose                 Verbose output
  --help                    Show this help

Examples:
  node test/config/test-runner.js --files backend/services/AtomPredictor.js
  node test/config/test-runner.js --layer unit
  node test/config/test-runner.js --all
      `);
      return;
    }

    if (args.includes('--integration-setup')) {
      await runner.setupIntegration();
      return;
    }

    if (args.includes('--all')) {
      const success = await runner.runAll();
      process.exit(success ? 0 : 1);
    }

    const layerIndex = args.indexOf('--layer');
    if (layerIndex !== -1 && args[layerIndex + 1]) {
      const layer = args[layerIndex + 1];
      if (Object.values(TEST_LAYERS).includes(layer)) {
        const success = await runner.runForLayer(layer);
        process.exit(success ? 0 : 1);
      } else {
        console.error(`❌ Invalid layer: ${layer}`);
        process.exit(1);
      }
    }

    const filesIndex = args.indexOf('--files');
    if (filesIndex !== -1 && args[filesIndex + 1]) {
      const files = args[filesIndex + 1].split(',');
      const success = await runner.runForFiles(files);
      process.exit(success ? 0 : 1);
    }

    // Default: run unit tests
    console.log('No specific options provided. Running unit tests...');
    const success = await runner.runForLayer(TEST_LAYERS.UNIT);
    process.exit(success ? 0 : 1);

  } catch (error) {
    console.error('❌ Test runner error:', error.message);
    runner.stopServer();
    process.exit(1);
  }
}

if (require.main === module) {
  main();
}

module.exports = TestRunner;