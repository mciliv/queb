#!/usr/bin/env node
/**
 * Output Capture Example
 * 
 * This example demonstrates how to use the output capture system
 * to automatically pass terminal output to AI commands for analysis.
 */

const OutputCapture = require('../src/core/OutputCapture');
const CommandRunner = require('../scripts/run-with-ai.js');

// ============================================================================
// EXAMPLE 1: Basic Output Capture
// ============================================================================

async function basicCaptureExample() {
  console.log('\n=== Example 1: Basic Output Capture ===\n');
  
  const capture = new OutputCapture({
    aiCommand: '~/Code/ai',
    autoTrigger: true,
    logger: console
  });
  
  // Simulate capturing output from a failing command
  console.log('Simulating a failing command...');
  
  const mockResult = {
    exitCode: 1,
    duration: 2500,
    output: [
      { type: 'stdout', content: 'Starting application...\n', timestamp: new Date() },
      { type: 'stdout', content: 'Loading configuration...\n', timestamp: new Date() }
    ],
    errors: [
      { type: 'stderr', content: 'Error: listen EADDRINUSE: address already in use :::8080\n', timestamp: new Date() },
      { type: 'stderr', content: '    at Server.setupListenHandle (net.js:1330:14)\n', timestamp: new Date() }
    ],
    logs: [],
    success: false
  };
  
  console.log('Mock result:', {
    exitCode: mockResult.exitCode,
    success: mockResult.success,
    outputLines: mockResult.output.length,
    errorLines: mockResult.errors.length
  });
  
  // Check if should trigger AI
  const shouldTrigger = capture.shouldTriggerAI(mockResult);
  console.log('Should trigger AI:', shouldTrigger);
  
  // Format for AI
  const formatted = capture.formatForAI(mockResult, {
    command: 'npm start',
    workingDirectory: '/Users/m/Code/queb',
    environment: 'development'
  });
  
  console.log('\nFormatted output for AI:');
  console.log('---');
  console.log(formatted.substring(0, 500) + '...');
  console.log('---');
}

// ============================================================================
// EXAMPLE 2: Command Runner Integration
// ============================================================================

async function commandRunnerExample() {
  console.log('\n=== Example 2: Command Runner Integration ===\n');
  
  const runner = new CommandRunner({
    aiCommand: '~/Code/ai',
    autoTrigger: true,
    forceAI: false,
    logger: console
  });
  
  // Test with a simple command that will fail
  console.log('Testing with a command that will fail...');
  
  try {
    const result = await runner.runWithAI('node', ['-e', 'console.log("Hello"); throw new Error("Test error");']);
    
    console.log('Command result:', {
      exitCode: result.command.exitCode,
      success: result.command.success,
      duration: result.command.duration,
      aiTriggered: result.ai !== null
    });
    
    if (result.ai) {
      console.log('AI analysis triggered:', result.ai.success);
    }
    
  } catch (error) {
    console.log('Command failed as expected:', error.message);
  }
}

// ============================================================================
// EXAMPLE 3: Real Command Testing
// ============================================================================

async function realCommandExample() {
  console.log('\n=== Example 3: Real Command Testing ===\n');
  
  const capture = new OutputCapture({
    aiCommand: 'echo', // Use echo instead of AI for testing
    autoTrigger: true,
    logger: console
  });
  
  // Test with a real command
  console.log('Testing with real command: ls /nonexistent');
  
  try {
    const result = await capture.captureCommand('ls /nonexistent');
    
    console.log('Command result:', {
      exitCode: result.exitCode,
      success: result.success,
      duration: result.duration,
      outputLines: result.output.length,
      errorLines: result.errors.length
    });
    
    // Show the actual output
    if (result.output.length > 0) {
      console.log('\nStandard output:');
      result.output.forEach(entry => {
        console.log(`[${entry.timestamp.toISOString()}] ${entry.content.trim()}`);
      });
    }
    
    if (result.errors.length > 0) {
      console.log('\nError output:');
      result.errors.forEach(entry => {
        console.log(`[${entry.timestamp.toISOString()}] ${entry.content.trim()}`);
      });
    }
    
    // Test AI triggering
    const shouldTrigger = capture.shouldTriggerAI(result);
    console.log('\nShould trigger AI:', shouldTrigger);
    
    if (shouldTrigger) {
      console.log('Would trigger AI analysis...');
      // In a real scenario, this would call the AI command
    }
    
  } catch (error) {
    console.log('Command execution error:', error.message);
  }
}

// ============================================================================
// EXAMPLE 4: Custom Error Patterns
// ============================================================================

async function customPatternsExample() {
  console.log('\n=== Example 4: Custom Error Patterns ===\n');
  
  const capture = new OutputCapture({
    aiCommand: '~/Code/ai',
    autoTrigger: true,
    triggerPatterns: [
      'custom error', 'special failure', 'unique issue'
    ],
    logger: console
  });
  
  // Test with custom error patterns
  const testCases = [
    {
      name: 'Standard error (should trigger)',
      output: 'Error: Something went wrong',
      errors: 'Failed to start'
    },
    {
      name: 'Custom error (should trigger)',
      output: 'custom error occurred',
      errors: ''
    },
    {
      name: 'Success message (should not trigger)',
      output: 'Application started successfully',
      errors: ''
    }
  ];
  
  for (const testCase of testCases) {
    const mockResult = {
      exitCode: 1,
      duration: 1000,
      output: [{ type: 'stdout', content: testCase.output, timestamp: new Date() }],
      errors: testCase.errors ? [{ type: 'stderr', content: testCase.errors, timestamp: new Date() }] : [],
      logs: [],
      success: false
    };
    
    const shouldTrigger = capture.shouldTriggerAI(mockResult);
    console.log(`${testCase.name}: ${shouldTrigger ? '✅ Would trigger AI' : '❌ Would not trigger AI'}`);
  }
}

// ============================================================================
// Run all examples
// ============================================================================

async function runAllExamples() {
  console.log('🚀 Output Capture Examples\n');
  console.log('This demonstrates how to capture terminal output and pass it to AI commands.\n');
  
  try {
    await basicCaptureExample();
    await commandRunnerExample();
    await realCommandExample();
    await customPatternsExample();
    
    console.log('\n✅ All examples completed successfully!\n');
    console.log('Key takeaways:');
    console.log('1. Output capture can monitor any command execution');
    console.log('2. AI analysis is triggered based on configurable patterns');
    console.log('3. Output is formatted for optimal AI consumption');
    console.log('4. The system is completely separate from error handling');
    console.log('5. Can be used with any AI command or script');
    console.log('6. Supports both real-time and batch processing');
    
  } catch (error) {
    console.error('❌ Error running examples:', error.message);
  }
  
  process.exit(0);
}

// Run if called directly
if (require.main === module) {
  runAllExamples();
}

module.exports = {
  basicCaptureExample,
  commandRunnerExample,
  realCommandExample,
  customPatternsExample
};