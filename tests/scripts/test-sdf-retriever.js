#!/usr/bin/env node

/**
 * SDF Retriever Test Runner
 * 
 * This script provides various ways to test the SDF retriever functionality:
 * - Unit tests for core logic
 * - Integration tests for end-to-end functionality  
 * - Performance tests for API calls
 * - Manual testing with specific compounds
 */

const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');

const CRYSTAL_SCRIPT = path.join(__dirname, '../../src/server/molecular-docking-research/crystal.py');
const TEST_DIR = path.join(__dirname, '..');

class SDFRetrieverTestRunner {
  constructor() {
    this.results = {
      unit: null,
      integration: null,
      manual: []
    };
  }

  async runCommand(command, args = [], options = {}) {
    return new Promise((resolve, reject) => {
      const child = spawn(command, args, {
        stdio: options.silent ? 'pipe' : 'inherit',
        cwd: options.cwd || TEST_DIR,
        timeout: options.timeout || 30000
      });

      let stdout = '';
      let stderr = '';

      if (options.silent) {
        child.stdout.on('data', (data) => stdout += data.toString());
        child.stderr.on('data', (data) => stderr += data.toString());
      }

      child.on('close', (code) => {
        resolve({ code, stdout: stdout.trim(), stderr: stderr.trim() });
      });

      child.on('error', reject);
    });
  }

  async runUnitTests() {
    console.log('🧪 Running SDF Retriever Unit Tests...');
    console.log('=' * 50);

    try {
      const result = await this.runCommand('npx', [
        'jest',
        'suites/unit/sdf-retriever.test.js',
        '--verbose'
      ]);

      this.results.unit = result;
      
      if (result.code === 0) {
        console.log('✅ Unit tests passed');
      } else {
        console.log('❌ Unit tests failed');
      }

      return result.code === 0;
    } catch (error) {
      console.error('💥 Unit test error:', error.message);
      return false;
    }
  }

  async runIntegrationTests() {
    console.log('\n🔗 Running SDF Retriever Integration Tests...');
    console.log('=' * 50);

    try {
      const result = await this.runCommand('npx', [
        'jest',
        'suites/integration/sdf-retriever.test.js',
        '--verbose',
        '--testTimeout=60000'
      ]);

      this.results.integration = result;
      
      if (result.code === 0) {
        console.log('✅ Integration tests passed');
      } else {
        console.log('❌ Integration tests failed');
      }

      return result.code === 0;
    } catch (error) {
      console.error('💥 Integration test error:', error.message);
      return false;
    }
  }

  async runManualTest(compound, description = '') {
    console.log(`\n🔍 Manual Test: ${compound} ${description ? '(' + description + ')' : ''}`);
    
    try {
      const tempDir = path.join(TEST_DIR, 'temp_manual_test');
      if (!fs.existsSync(tempDir)) {
        fs.mkdirSync(tempDir, { recursive: true });
      }

      const result = await this.runCommand('python', [
        CRYSTAL_SCRIPT,
        compound,
        '--dir',
        tempDir
      ], { silent: true, timeout: 45000 });

      const testResult = {
        compound,
        description,
        success: result.code === 0,
        output: result.stdout,
        error: result.stderr,
        timestamp: new Date().toISOString()
      };

      this.results.manual.push(testResult);

      if (result.code === 0) {
        console.log(`  ✅ SUCCESS: ${result.stdout.split('\n')[0]}`);
        
        // Check solution type
        if (result.stdout.includes('Local file system')) {
          console.log('  📁 Source: Local file system');
        } else if (result.stdout.includes('API download')) {
          console.log('  🌐 Source: API download');
        }
        
        // Show file info
        const sizeMatch = result.stdout.match(/File size: (\d+) bytes/);
        if (sizeMatch) {
          console.log(`  📄 Size: ${sizeMatch[1]} bytes`);
        }
      } else {
        console.log(`  ❌ FAILED: ${result.stderr || 'No SDF found'}`);
      }

      // Cleanup
      try {
        fs.rmSync(tempDir, { recursive: true, force: true });
      } catch (e) {
        // Ignore cleanup errors
      }

      return testResult;
    } catch (error) {
      console.log(`  💥 ERROR: ${error.message}`);
      return {
        compound,
        description,
        success: false,
        error: error.message,
        timestamp: new Date().toISOString()
      };
    }
  }

  async runComprehensiveTestSuite() {
    console.log('\n🎯 Running Comprehensive SDF Retriever Test Suite...');
    console.log('=' * 50);

    try {
      const result = await this.runCommand('python', [
        CRYSTAL_SCRIPT,
        '--test-suite'
      ], { timeout: 120000 });

      if (result.code === 0) {
        console.log('✅ Comprehensive test suite completed successfully');
        
        // Parse results
        const successMatch = result.stdout.match(/Overall Success Rate: ([\d.]+)%/);
        if (successMatch) {
          console.log(`📈 Success Rate: ${successMatch[1]}%`);
        }

        const localMatch = result.stdout.match(/Local file retrieval: (\d+) cases/);
        const apiMatch = result.stdout.match(/API-based retrieval: (\d+) cases/);
        
        if (localMatch) console.log(`📁 Local retrievals: ${localMatch[1]}`);
        if (apiMatch) console.log(`🌐 API retrievals: ${apiMatch[1]}`);
      } else {
        console.log('❌ Comprehensive test suite failed');
      }

      return result.code === 0;
    } catch (error) {
      console.error('💥 Comprehensive test error:', error.message);
      return false;
    }
  }

  async runPerformanceTest() {
    console.log('\n⚡ Running Performance Tests...');
    console.log('=' * 40);

    const testCases = [
      { compound: 'CCO', expected: 'local', description: 'Local file (fast)' },
      { compound: 'caffeine', expected: 'api', description: 'API call (slower)' },
      { compound: 'aspirin', expected: 'api', description: 'API call (slower)' }
    ];

    const results = [];

    for (const testCase of testCases) {
      const startTime = Date.now();
      const result = await this.runManualTest(testCase.compound, testCase.description);
      const endTime = Date.now();
      const duration = endTime - startTime;

      results.push({
        ...result,
        duration,
        expected: testCase.expected
      });

      console.log(`  ⏱️  Duration: ${duration}ms`);
    }

    // Performance summary
    console.log('\n📊 Performance Summary:');
    const avgLocal = results.filter(r => r.output && r.output.includes('Local file system'))
                           .reduce((sum, r) => sum + r.duration, 0) / results.filter(r => r.output && r.output.includes('Local file system')).length || 0;
    const avgAPI = results.filter(r => r.output && r.output.includes('API download'))
                         .reduce((sum, r) => sum + r.duration, 0) / results.filter(r => r.output && r.output.includes('API download')).length || 0;

    if (avgLocal > 0) console.log(`  📁 Average local retrieval: ${Math.round(avgLocal)}ms`);
    if (avgAPI > 0) console.log(`  🌐 Average API retrieval: ${Math.round(avgAPI)}ms`);

    return results;
  }

  generateReport() {
    console.log('\n📋 Test Report Summary');
    console.log('=' * 50);

    if (this.results.unit) {
      console.log(`Unit Tests: ${this.results.unit.code === 0 ? '✅ PASSED' : '❌ FAILED'}`);
    }

    if (this.results.integration) {
      console.log(`Integration Tests: ${this.results.integration.code === 0 ? '✅ PASSED' : '❌ FAILED'}`);
    }

    if (this.results.manual.length > 0) {
      const passed = this.results.manual.filter(r => r.success).length;
      const total = this.results.manual.length;
      console.log(`Manual Tests: ${passed}/${total} passed (${Math.round(passed/total*100)}%)`);
      
      console.log('\nManual Test Details:');
      this.results.manual.forEach(result => {
        console.log(`  ${result.success ? '✅' : '❌'} ${result.compound} ${result.description ? '(' + result.description + ')' : ''}`);
      });
    }

    console.log(`\nReport generated at: ${new Date().toISOString()}`);
  }
}

// CLI Interface
async function main() {
  const args = process.argv.slice(2);
  const runner = new SDFRetrieverTestRunner();

  if (args.length === 0) {
    console.log('SDF Retriever Test Runner');
    console.log('Usage:');
    console.log('  node test-sdf-retriever.js [command]');
    console.log('');
    console.log('Commands:');
    console.log('  unit          Run unit tests');
    console.log('  integration   Run integration tests');
    console.log('  manual        Run manual tests with common compounds');
    console.log('  comprehensive Run the built-in comprehensive test suite');
    console.log('  performance   Run performance benchmarks');
    console.log('  all           Run all tests');
    console.log('  test <compound>  Test specific compound');
    console.log('');
    console.log('Examples:');
    console.log('  node test-sdf-retriever.js all');
    console.log('  node test-sdf-retriever.js test caffeine');
    console.log('  node test-sdf-retriever.js performance');
    return;
  }

  const command = args[0];

  try {
    switch (command) {
      case 'unit':
        await runner.runUnitTests();
        break;

      case 'integration':
        await runner.runIntegrationTests();
        break;

      case 'manual':
        console.log('Running manual tests with common compounds...');
        await runner.runManualTest('CCO', 'Ethanol - should be local');
        await runner.runManualTest('caffeine', 'Should use API');
        await runner.runManualTest('aspirin', 'Should use API');
        await runner.runManualTest('XYZ123', 'Invalid - should fail');
        break;

      case 'comprehensive':
        await runner.runComprehensiveTestSuite();
        break;

      case 'performance':
        await runner.runPerformanceTest();
        break;

      case 'all':
        await runner.runUnitTests();
        await runner.runIntegrationTests();
        await runner.runComprehensiveTestSuite();
        await runner.runPerformanceTest();
        break;

      case 'test':
        if (args[1]) {
          await runner.runManualTest(args[1], 'Manual test');
        } else {
          console.log('Please specify a compound to test');
        }
        break;

      default:
        console.log(`Unknown command: ${command}`);
        console.log('Use "node test-sdf-retriever.js" for help');
        process.exit(1);
    }

    runner.generateReport();
  } catch (error) {
    console.error('Test runner error:', error);
    process.exit(1);
  }
}

if (require.main === module) {
  main();
}

module.exports = SDFRetrieverTestRunner;
