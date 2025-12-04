#!/usr/bin/env node
/**
 * External Script Integration Example
 * 
 * This example shows how to use the external script trigger system
 * for automatic error resolution. It demonstrates:
 * 
 * 1. How to configure external scripts
 * 2. How to test the system
 * 3. How to create custom error handlers
 * 4. How to integrate with existing error handling
 */

const ExternalScriptTrigger = require('../src/core/ExternalScriptTrigger');
const ErrorHandler = require('../src/core/ErrorHandler');

// ============================================================================
// EXAMPLE 1: Basic Configuration and Testing
// ============================================================================

async function basicExample() {
  console.log('\n=== Example 1: Basic External Script Trigger ===\n');
  
  // Create a script trigger with custom configuration
  const scriptTrigger = new ExternalScriptTrigger({
    enabled: true,
    scripts: [
      {
        name: 'ai-script',
        path: '~/Code/ai',
        enabled: true,
        priority: 1,
        timeout: 30000,
        retries: 2
      },
      {
        name: 'echo-test',
        path: 'echo',
        enabled: true,
        priority: 2,
        timeout: 5000,
        retries: 1,
        args: ['Error occurred:', '${ERROR_MESSAGE}']
      }
    ],
    minSeverity: 'medium',
    logger: console
  });
  
  // Simulate an error
  const mockError = {
    message: 'listen EADDRINUSE: address already in use :::8080',
    code: 'EADDRINUSE',
    severity: 'medium',
    category: 'network',
    timestamp: new Date().toISOString(),
    stack: 'Error: listen EADDRINUSE: address already in use :::8080\n    at Server.setupListenHandle'
  };
  
  console.log('Simulating error:', mockError.message);
  console.log('Should trigger scripts:', scriptTrigger.shouldTriggerScripts(mockError));
  
  // Execute scripts (this would normally be called by ErrorHandler)
  const result = await scriptTrigger.executeScripts(mockError);
  console.log('Script execution result:', result);
}

// ============================================================================
// EXAMPLE 2: Integration with ErrorHandler
// ============================================================================

async function errorHandlerExample() {
  console.log('\n=== Example 2: ErrorHandler Integration ===\n');
  
  // Create error handler with external script integration
  const errorHandler = new ErrorHandler();
  errorHandler.initialize(console);
  
  // Simulate different types of errors
  const errors = [
    {
      error: new Error('listen EADDRINUSE: address already in use :::8080'),
      context: { category: 'network', userFacing: true }
    },
    {
      error: new Error('Cannot read property "length" of undefined'),
      context: { category: 'application', userFacing: false }
    },
    {
      error: new Error('Invalid user input format'),
      context: { category: 'validation', userFacing: true }
    }
  ];
  
  for (const { error, context } of errors) {
    console.log(`\nHandling error: ${error.message}`);
    const result = await errorHandler.handle(error, context);
    
    console.log('Error handled:', {
      message: result.message,
      code: result.code,
      recoverable: result.recoverable,
      externalScripts: result.externalScripts ? 'Triggered' : 'Not triggered'
    });
    
    if (result.externalScripts) {
      console.log('Script results:', result.externalScripts.scripts.map(s => ({
        name: s.script,
        success: s.success,
        duration: s.duration
      })));
    }
  }
}

// ============================================================================
// EXAMPLE 3: Custom Script Creation
// ============================================================================

async function customScriptExample() {
  console.log('\n=== Example 3: Custom Script Creation ===\n');
  
  // Create a simple error handler script
  const customScript = `#!/bin/bash
# Custom Error Handler Script
# This script receives error context via environment variables

echo "=== Custom Error Handler ==="
echo "Error: $ERROR_MESSAGE"
echo "Code: $ERROR_CODE"
echo "Severity: $ERROR_SEVERITY"
echo "Category: $ERROR_CATEGORY"
echo "Project: $PROJECT_NAME ($PROJECT_VERSION)"
echo "Platform: $PLATFORM ($ARCH)"
echo "Environment: Development=$IS_DEVELOPMENT, Test=$IS_TEST, Production=$IS_PRODUCTION"
echo "Timestamp: $ERROR_TIMESTAMP"
echo ""

# Simple error-specific handling
case "$ERROR_CODE" in
  "EADDRINUSE")
    echo "🔧 Suggested fix: Kill process using port and restart"
    echo "   lsof -ti:8080 | xargs kill -9"
    ;;
  "ENOTFOUND")
    echo "🔧 Suggested fix: Check network connection and DNS"
    ;;
  "MODULE_NOT_FOUND")
    echo "🔧 Suggested fix: Run npm install or check module path"
    ;;
  *)
    echo "🔧 General suggestion: Check logs and retry"
    ;;
esac

echo "=== End Custom Error Handler ==="
`;
  
  const fs = require('fs');
  const path = require('path');
  
  // Write the script to a temporary location
  const scriptPath = path.join(__dirname, 'custom-error-handler.sh');
  fs.writeFileSync(scriptPath, customScript);
  fs.chmodSync(scriptPath, '755');
  
  console.log('Created custom script at:', scriptPath);
  
  // Create script trigger with custom script
  const scriptTrigger = new ExternalScriptTrigger({
    enabled: true,
    scripts: [
      {
        name: 'custom-handler',
        path: scriptPath,
        enabled: true,
        priority: 1,
        timeout: 10000,
        retries: 1
      }
    ],
    logger: console
  });
  
  // Test with an error
  const mockError = {
    message: 'listen EADDRINUSE: address already in use :::8080',
    code: 'EADDRINUSE',
    severity: 'medium',
    category: 'network',
    timestamp: new Date().toISOString()
  };
  
  console.log('Testing custom script with EADDRINUSE error...');
  const result = await scriptTrigger.executeScripts(mockError);
  
  if (result && result.scripts.length > 0) {
    console.log('Custom script output:');
    console.log(result.scripts[0].output);
  }
  
  // Cleanup
  fs.unlinkSync(scriptPath);
  console.log('Cleaned up custom script');
}

// ============================================================================
// EXAMPLE 4: Configuration Management
// ============================================================================

async function configurationExample() {
  console.log('\n=== Example 4: Configuration Management ===\n');
  
  const scriptTrigger = new ExternalScriptTrigger({
    enabled: true,
    scripts: [
      { name: 'script1', path: 'echo', enabled: true, priority: 1 },
      { name: 'script2', path: 'echo', enabled: false, priority: 2 }
    ],
    logger: console
  });
  
  console.log('Initial configuration:');
  console.log(JSON.stringify(scriptTrigger.getConfig(), null, 2));
  
  // Add a new script
  scriptTrigger.addScript({
    name: 'new-script',
    path: 'echo',
    enabled: true,
    priority: 1.5,
    args: ['New script added!']
  });
  
  console.log('\nAfter adding new script:');
  console.log(JSON.stringify(scriptTrigger.getConfig(), null, 2));
  
  // Enable/disable scripts
  scriptTrigger.setScriptEnabled('script2', true);
  scriptTrigger.setScriptEnabled('script1', false);
  
  console.log('\nAfter enabling/disabling scripts:');
  console.log(JSON.stringify(scriptTrigger.getConfig(), null, 2));
}

// ============================================================================
// Run all examples
// ============================================================================

async function runAllExamples() {
  console.log('🚀 External Script Integration Examples\n');
  console.log('This demonstrates how to use external scripts for automatic error resolution.\n');
  
  try {
    await basicExample();
    await errorHandlerExample();
    await customScriptExample();
    await configurationExample();
    
    console.log('\n✅ All examples completed successfully!\n');
    console.log('Key takeaways:');
    console.log('1. External scripts are triggered automatically on errors');
    console.log('2. Scripts receive error context via environment variables');
    console.log('3. Multiple scripts can be configured with priorities');
    console.log('4. Rate limiting prevents script spam');
    console.log('5. Scripts can be enabled/disabled dynamically');
    console.log('6. The system is completely abstracted from specific tools');
    
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
  basicExample,
  errorHandlerExample,
  customScriptExample,
  configurationExample
};