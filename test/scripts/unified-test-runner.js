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

  async isFrontendAvailable(url = 'http://localhost:3001') {
    try {
      const response = await fetch(url, { method: 'HEAD' });
      // Consider any non-network response as available (even 404)
      return response && typeof response.status === 'number';
    } catch (_) {
      return false;
    }
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
    try {
      const success = await this.runCommand('npx', [
        'jest', 
        '--testPathPattern=unit.test.js',
        '--silent',
        '--reporters=default'
      ]);
      this.results.unit = success;
      // Unit test results recorded
      return success;
    } catch (error) {
      this.results.unit = false;
      return false;
    }
  }

  async runVisualTests() {
    try {
      const available = await this.isFrontendAvailable('http://localhost:3001');
      if (!available) {
        this.results.visual = 'skipped';
        return true;
      }
      // Use AutoTabConnector for seamless tab management
      const tabInfo = await AutoTabConnector.getOrCreateTab();
      
      const success = await this.runCommand('npx', [
        'jest',
        'test/suites/integration/visual-interface.test.js',
        '--testTimeout=30000',
        '--silent'
      ]);
      
      this.results.visual = success;
      // Visual test results recorded
      return success;
    } catch (error) {
      this.results.visual = false;
      return false;
    }
  }

  async runMolecularTests(customCases = null) {
    const testCases = customCases || ['water', 'ethanol', 'coffee'];
    
    try {
      const available = await this.isFrontendAvailable('http://localhost:3001');
      if (!available) {
        this.results.molecular = 'skipped';
        return true;
      }
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
    try {
      if (!process.env.OPENAI_API_KEY) {
        this.results.pipeline = 'skipped';
        return true;
      }

      const success = await this.runCommand('npx', [
        'jest',
        'test/suites/integration/full-pipeline-visualization.test.js',
        '--testTimeout=60000',
        '--silent'
      ]);
      
      this.results.pipeline = success;
      // Pipeline test results recorded
      return success;
    } catch (error) {
      this.results.pipeline = false;
      return false;
    }
  }

  async runPersistenceTests() {
    try {
      const available = await this.isFrontendAvailable('http://localhost:3001');
      if (!available) {
        this.results.persistence = 'skipped';
        return true;
      }
      const success = await this.runCommand('npx', [
        'jest',
        'test/suites/integration/persistent-tab-tests.test.js',
        '--testTimeout=30000',
        '--silent'
      ]);
      
      this.results.persistence = success;
      // Persistence test results recorded
      return success;
    } catch (error) {
      this.results.persistence = false;
      return false;
    }
  }

  async runApiTests() {
    try {
      // Check if OpenAI API key is available
      if (!process.env.OPENAI_API_KEY) {
        this.results.integration = 'skipped';
        return true;
      }

      const success = await this.runCommand('npx', [
        'jest',
        'test/suites/integration/molecular-accuracy.test.js',
        '--testTimeout=60000',
        '--silent'
      ]);
      
      this.results.integration = success;
      // API test results recorded
      return success;
    } catch (error) {
      this.results.integration = false;
      return false;
    }
  }

  async runCommandLineTests() {
    try {
      const success = await this.runCommand('node', [
        'test/scripts/test-molecular-ui.js'
      ], { silent: true });
      
      return success;
    } catch (error) {
      return false;
    }
  }

  async runDemoTests() {
    try {
      const success = await this.runCommand('node', [
        'test/scripts/test-visual-demo.js'
      ], { silent: true });
      
      return success;
    } catch (error) {
      return false;
    }
  }

  async runFakeDataTests() {
    try {
      const success = await this.runCommand('node', [
        'test/scripts/test-inject-fake-data.js'
      ], { silent: true });
      
      return success;
    } catch (error) {
      return false;
    }
  }

  async runInjectionTests() {
    try {
      const success = await this.runCommand('npx', [
        'jest',
        'test/suites/integration/auto-inject-tests.test.js',
        '--testTimeout=60000',
        '--silent'
      ]);
      
      return success;
    } catch (error) {
      return false;
    }
  }

  printSummary() {
    const considered = Object.values(this.results).filter(r => r !== null && r !== 'skipped');
    const totalTests = considered.length;
    const passedTests = considered.filter(r => r === true).length;
    const failedTests = considered.filter(r => r === false).length;
    
    // Always show summary
    console.log(`\n${passedTests}/${totalTests} test suites passed`);
    
    if (failedTests > 0) {
      console.log('\n❌ Failed:');
      if (this.results.unit === false) console.log('  Unit Tests');
      if (this.results.visual === false) console.log('  Visual Interface Tests');
      if (this.results.molecular === false) console.log('  Molecular Visualization Tests');
      if (this.results.pipeline === false) console.log('  Full Pipeline Tests');
      if (this.results.persistence === false) console.log('  Tab Persistence Tests');
      if (this.results.integration === false) console.log('  API Integration Tests');
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
  runner.run().catch(error => {
    console.error('❌ Test runner error:', error.message);
    process.exit(1);
  });
}

module.exports = UnifiedTestRunner;