#!/usr/bin/env node

const { spawn } = require('child_process');
const AutoTabConnector = require('../utils/auto-tab-connector');

class UnifiedTestRunner {
  constructor() {
    this.results = {
      unit: null,
      integration: null,
      visual: null,
      molecular: null
    };
  }

  async runCommand(command, args = [], options = {}) {
    return new Promise((resolve, reject) => {
      const child = spawn(command, args, {
        stdio: options.silent ? 'ignore' : 'inherit',
        shell: true,
        ...options
      });

      child.on('close', (code) => {
        resolve(code === 0);
      });

      child.on('error', (error) => {
        reject(error);
      });
    });
  }

  async runUnitTests() {
    console.log('🧪 Running unit tests...');
    try {
      const success = await this.runCommand('npx', [
        'jest', 
        '--testPathPattern=unit.test.js',
        '--verbose',
        '--silent'
      ]);
      this.results.unit = success;
      if (success) {
        console.log('✅ Unit tests passed');
      } else {
        console.log('❌ Unit tests failed');
      }
      return success;
    } catch (error) {
      console.log('❌ Unit tests error:', error.message);
      this.results.unit = false;
      return false;
    }
  }

  async runVisualTests() {
    console.log('🖥️  Running visual interface tests...');
    try {
      // Use AutoTabConnector for seamless tab management
      const tabInfo = await AutoTabConnector.getOrCreateTab();
      console.log(tabInfo.reused ? '♻️  Reusing existing Chrome tab' : '🌟 Created new Chrome tab');
      
      const success = await this.runCommand('npx', [
        'jest',
        'test/integration/visual-interface.test.js',
        '--testTimeout=30000',
        '--silent'
      ]);
      
      this.results.visual = success;
      if (success) {
        console.log('✅ Visual tests passed');
      } else {
        console.log('❌ Visual tests failed');
      }
      return success;
    } catch (error) {
      console.log('❌ Visual tests error:', error.message);
      this.results.visual = false;
      return false;
    }
  }

  async runMolecularTests() {
    console.log('🧬 Running molecular visualization tests...');
    
    const testCases = ['water', 'ethanol', 'coffee'];
    
    try {
      for (const testCase of testCases) {
        console.log(`  🔬 Testing: ${testCase}`);
        await AutoTabConnector.executeTestSilently(testCase);
        await new Promise(resolve => setTimeout(resolve, 3000)); // Brief wait
      }
      
      this.results.molecular = true;
      console.log('✅ Molecular tests completed');
      return true;
    } catch (error) {
      console.log('❌ Molecular tests error:', error.message);
      this.results.molecular = false;
      return false;
    }
  }

  async runApiTests() {
    console.log('🔗 Running API integration tests...');
    try {
      // Check if OpenAI API key is available
      if (!process.env.OPENAI_API_KEY) {
        console.log('⚠️  Skipping API tests - no OPENAI_API_KEY');
        this.results.integration = 'skipped';
        return true;
      }

      const success = await this.runCommand('npx', [
        'jest',
        'test/integration/molecular-accuracy.test.js',
        '--testTimeout=60000',
        '--silent'
      ]);
      
      this.results.integration = success;
      if (success) {
        console.log('✅ API integration tests passed');
      } else {
        console.log('❌ API integration tests failed');
      }
      return success;
    } catch (error) {
      console.log('❌ API tests error:', error.message);
      this.results.integration = false;
      return false;
    }
  }

  printSummary() {
    console.log('\n📊 Test Summary');
    console.log('================');
    
    const statusIcon = (result) => {
      if (result === true) return '✅';
      if (result === 'skipped') return '⏭️ ';
      if (result === false) return '❌';
      return '❓';
    };

    console.log(`${statusIcon(this.results.unit)} Unit Tests`);
    console.log(`${statusIcon(this.results.visual)} Visual Interface Tests`);
    console.log(`${statusIcon(this.results.molecular)} Molecular Visualization Tests`);
    console.log(`${statusIcon(this.results.integration)} API Integration Tests`);

    const failedTests = Object.values(this.results).filter(r => r === false).length;
    const totalTests = Object.values(this.results).filter(r => r !== null).length;
    
    if (failedTests === 0) {
      console.log('\n🎉 All tests passed!');
      console.log('💡 Chrome tab remains open for continued testing');
    } else {
      console.log(`\n⚠️  ${failedTests}/${totalTests} test suites failed`);
    }
  }

  async run() {
    console.log('🚀 Unified Test Runner');
    console.log('======================');
    console.log('Running all molecular analysis tests...\n');

    // Run tests in parallel where possible
    const unitSuccess = await this.runUnitTests();
    
    // Only continue with integration tests if unit tests pass
    if (unitSuccess) {
      // Run visual and molecular tests in parallel
      await Promise.all([
        this.runVisualTests(),
        this.runMolecularTests()
      ]);
      
      // Run API tests last
      await this.runApiTests();
    } else {
      console.log('❌ Unit tests failed - skipping integration tests');
    }

    this.printSummary();
    
    // Exit code based on critical test failures
    const criticalFailures = [this.results.unit, this.results.visual].filter(r => r === false).length;
    process.exit(criticalFailures > 0 ? 1 : 0);
  }
}

// Run if called directly
if (require.main === module) {
  const runner = new UnifiedTestRunner();
  runner.run().catch(console.error);
}

module.exports = UnifiedTestRunner;