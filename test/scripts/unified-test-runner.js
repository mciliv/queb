#!/usr/bin/env node

const { spawn } = require('child_process');
const AutoTabConnector = require('../utils/auto-tab-connector');

class UnifiedTestRunner {
  constructor() {
    this.results = {
      unit: null,
      integration: null,
      visual: null,
      molecular: null,
      pipeline: null,
      persistence: null
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

  async runMolecularTests(customCases = null) {
    const testCases = customCases || ['water', 'ethanol', 'coffee'];
    
    try {
      for (const testCase of testCases) {
        const page = await AutoTabConnector.executeTestSilently(testCase);
        await new Promise(resolve => setTimeout(resolve, 6000));
        
        const screenshotPath = `test/screenshots/test_${testCase.replace(/\s+/g, '_')}.png`;
        await page.screenshot({ path: screenshotPath, fullPage: true });
      }
      
      this.results.molecular = true;
      return true;
    } catch (error) {
      this.results.molecular = false;
      return false;
    }
  }

  async runPipelineTests() {
    console.log('🧬 Running full pipeline tests...');
    try {
      if (!process.env.OPENAI_API_KEY) {
        console.log('⚠️  Skipping pipeline tests - no OPENAI_API_KEY');
        this.results.pipeline = 'skipped';
        return true;
      }

      const success = await this.runCommand('npx', [
        'jest',
        'test/integration/full-pipeline-visualization.test.js',
        '--testTimeout=60000',
        '--silent'
      ]);
      
      this.results.pipeline = success;
      if (success) {
        console.log('✅ Pipeline tests passed');
      } else {
        console.log('❌ Pipeline tests failed');
      }
      return success;
    } catch (error) {
      console.log('❌ Pipeline tests error:', error.message);
      this.results.pipeline = false;
      return false;
    }
  }

  async runPersistenceTests() {
    console.log('🔄 Running persistence/tab management tests...');
    try {
      const success = await this.runCommand('npx', [
        'jest',
        'test/integration/persistent-tab-tests.test.js',
        '--testTimeout=30000',
        '--silent'
      ]);
      
      this.results.persistence = success;
      if (success) {
        console.log('✅ Persistence tests passed');
      } else {
        console.log('❌ Persistence tests failed');
      }
      return success;
    } catch (error) {
      console.log('❌ Persistence tests error:', error.message);
      this.results.persistence = false;
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

  async runCommandLineTests() {
    console.log('💻 Running command line molecular tests...');
    try {
      const success = await this.runCommand('node', [
        'test/scripts/test-molecular-ui.js'
      ], { silent: true });
      
      return success;
    } catch (error) {
      console.log('❌ Command line tests error:', error.message);
      return false;
    }
  }

  async runDemoTests() {
    console.log('🎭 Running visual demo tests...');
    try {
      const success = await this.runCommand('node', [
        'test/scripts/test-visual-demo.js'
      ], { silent: true });
      
      return success;
    } catch (error) {
      console.log('❌ Demo tests error:', error.message);
      return false;
    }
  }

  async runFakeDataTests() {
    console.log('🎯 Running fake data injection tests...');
    try {
      const success = await this.runCommand('node', [
        'test/scripts/test-inject-fake-data.js'
      ], { silent: true });
      
      return success;
    } catch (error) {
      console.log('❌ Fake data tests error:', error.message);
      return false;
    }
  }

  async runInjectionTests() {
    console.log('💉 Running visual injection tests...');
    try {
      const success = await this.runCommand('npx', [
        'jest',
        'test/integration/auto-inject-tests.test.js',
        '--testTimeout=60000',
        '--silent'
      ]);
      
      return success;
    } catch (error) {
      console.log('❌ Injection tests error:', error.message);
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
    console.log(`${statusIcon(this.results.pipeline)} Full Pipeline Tests`);
    console.log(`${statusIcon(this.results.persistence)} Tab Persistence Tests`);
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
    const args = process.argv.slice(2);
    
    // Handle specific test types
    if (args.includes('--visual')) {
      await this.runVisualTests();
      await this.runInjectionTests();
      this.printSummary();
      return;
    }
    
    if (args.includes('--molecular')) {
      await this.runMolecularTests();
      this.printSummary();
      return;
    }
    
    if (args.includes('--pipeline')) {
      await this.runPipelineTests();
      this.printSummary();
      return;
    }
    
    if (args.includes('--persistence')) {
      await this.runPersistenceTests();
      this.printSummary();
      return;
    }
    
    if (args.includes('--demo')) {
      await this.runDemoTests();
      this.printSummary();
      return;
    }
    
    if (args.includes('--fake-data')) {
      await this.runFakeDataTests();
      this.printSummary();
      return;
    }
    
    if (args.includes('--cli')) {
      await this.runCommandLineTests();
      this.printSummary();
      return;
    }
    
    if (args.includes('--quick')) {
      const testCases = args.filter(arg => !arg.startsWith('--'));
      await this.runMolecularTests(testCases.length > 0 ? testCases : ['water', 'ethanol']);
      this.printSummary();
      return;
    }
    
    // Quick mode: just run molecular tests with custom cases (legacy support)
    if (args.length > 0 && !args[0].startsWith('--')) {
      for (const testCase of args) {
        const page = await AutoTabConnector.executeTestSilently(testCase);
        await new Promise(resolve => setTimeout(resolve, 6000));
        const screenshotPath = `test/screenshots/test_${testCase.replace(/\s+/g, '_')}.png`;
        await page.screenshot({ path: screenshotPath, fullPage: true });
      }
      return;
    }

    // Full test suite
    console.log('🚀 Running comprehensive test suite...\n');
    
    const unitSuccess = await this.runUnitTests();
    
    if (unitSuccess) {
      // Run visual and molecular tests in parallel
      await Promise.all([
        this.runVisualTests(),
        this.runMolecularTests(),
        this.runPersistenceTests()
      ]);
      
      // Run API-dependent tests
      await Promise.all([
        this.runPipelineTests(),
        this.runApiTests()
      ]);
    }

    this.printSummary();
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