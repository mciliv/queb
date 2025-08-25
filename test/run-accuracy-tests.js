#!/usr/bin/env node
// test/run-accuracy-tests.js - Test runner for molecular accuracy validation
// Usage: node test/run-accuracy-tests.js [--api-key=YOUR_KEY] [--material=specific_material]

const { spawn } = require('child_process');
const path = require('path');

// Parse command line arguments
const args = process.argv.slice(2);
const options = {
  apiKey: null,
  material: null,
  verbose: false,
  suite: 'all' // all, unit, integration
};

args.forEach(arg => {
  if (arg.startsWith('--api-key=')) {
    options.apiKey = arg.split('=')[1];
  } else if (arg.startsWith('--material=')) {
    options.material = arg.split('=')[1];
  } else if (arg === '--verbose') {
    options.verbose = true;
  } else if (arg.startsWith('--suite=')) {
    options.suite = arg.split('=')[1];
  } else if (arg === '--help') {
    showHelp();
    process.exit(0);
  }
});

function showHelp() {
  console.log(`
Molecular Accuracy Test Runner

Usage: node test/run-accuracy-tests.js [options]

Options:
  --api-key=KEY        OpenAI API key for integration tests
  --material=NAME      Test specific material only
  --suite=TYPE         Test suite: all, unit, integration (default: all)
  --verbose            Show detailed output
  --help               Show this help

Examples:
  node test/run-accuracy-tests.js
  node test/run-accuracy-tests.js --api-key=sk-xxx --material=water
  node test/run-accuracy-tests.js --suite=unit
  node test/run-accuracy-tests.js --suite=integration --api-key=sk-xxx

Environment Variables:
  OPENAI_API_KEY       OpenAI API key (alternative to --api-key)
`);
}

async function runTests() {
  console.log('🧪 Starting Molecular Accuracy Tests...\n');

  // Set up environment
  if (options.apiKey) {
    process.env.OPENAI_API_KEY = options.apiKey;
  }

  if (!process.env.OPENAI_API_KEY && (options.suite === 'all' || options.suite === 'integration')) {
    console.log('⚠️  Warning: No OpenAI API key set. Integration tests will be skipped.');
    console.log('   Set OPENAI_API_KEY environment variable or use --api-key option\n');
  }

  // Determine which tests to run
  const testFiles = [];
  
  if (options.suite === 'all' || options.suite === 'unit') {
    testFiles.push('test/unit/prompt-accuracy.test.js');
  }
  
  if (options.suite === 'all' || options.suite === 'integration') {
    testFiles.push('test/integration/molecular-accuracy.test.js');
  }

  // Run each test suite
  for (const testFile of testFiles) {
    console.log(`\n📋 Running ${testFile}...`);
    
    const jestArgs = [
      '--config', 'test/jest.config.js',
      testFile,
      '--verbose',
      '--colors',
      '--no-cache'
    ];

    // Add material filter if specified
    if (options.material) {
      jestArgs.push('--testNamePattern', options.material);
    }

    try {
      await runJest(jestArgs);
      console.log(`✅ ${testFile} completed successfully`);
    } catch (error) {
      console.error(`❌ ${testFile} failed:`, error.message);
      if (!options.verbose) {
        console.log('   Use --verbose for detailed output');
      }
    }
  }

  console.log('\n🎯 Molecular accuracy testing completed!');
  console.log('\n📊 To view detailed results, check the test output above.');
  console.log('💡 Tip: Use --material=water to test specific materials');
}

function runJest(args) {
  return new Promise((resolve, reject) => {
    const jest = spawn('npx', ['jest', ...args], {
      stdio: options.verbose ? 'inherit' : 'pipe',
      cwd: path.resolve(__dirname, '..')
    });

    let output = '';
    
    if (!options.verbose) {
      jest.stdout?.on('data', (data) => {
        output += data.toString();
      });
      
      jest.stderr?.on('data', (data) => {
        output += data.toString();
      });
    }

    jest.on('close', (code) => {
      if (code === 0) {
        resolve();
      } else {
        if (!options.verbose && output) {
          console.log(output);
        }
        reject(new Error(`Jest exited with code ${code}`));
      }
    });

    jest.on('error', (error) => {
      reject(error);
    });
  });
}

// Run if called directly
if (require.main === module) {
  runTests().catch(error => {
    console.error('Test runner failed:', error.message);
    process.exit(1);
  });
}

module.exports = { runTests, options }; 