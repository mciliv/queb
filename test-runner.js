#!/usr/bin/env node

const { spawn } = require('child_process');
const path = require('path');

function runTest(testFile) {
  console.log(`🧪 Running tests...`);

  if (testFile) {
    console.log(`Running specific test: ${testFile}`);
    // For now, just run the test file directly if it exists
    const testPath = path.resolve(testFile);
    try {
      require(testPath);
    } catch (error) {
      console.error(`Failed to run test ${testFile}:`, error.message);
      process.exit(1);
    }
    return;
  }

  // Run all test scripts
  const testScripts = [
    'test-env-loading.js',
    'test-unified-ai-service.js',
    'test-completion-api.js'
  ];

  console.log('Running all test scripts...');

  let completed = 0;
  const total = testScripts.length;

  function runNext() {
    if (completed >= total) {
      console.log(`✅ All ${total} test scripts completed!`);
      return;
    }

    const script = testScripts[completed];
    console.log(`\n📋 Running ${script} (${completed + 1}/${total})`);

    const child = spawn('node', [script], {
      stdio: 'inherit',
      cwd: process.cwd()
    });

    child.on('close', (code) => {
      if (code === 0) {
        console.log(`✅ ${script} passed`);
      } else {
        console.log(`❌ ${script} failed with code ${code}`);
      }
      completed++;
      runNext();
    });

    child.on('error', (error) => {
      console.error(`❌ Failed to run ${script}:`, error.message);
      completed++;
      runNext();
    });
  }

  runNext();
}

module.exports = runTest;

// If called directly
if (require.main === module) {
  const testFile = process.argv[2];
  runTest(testFile);
}