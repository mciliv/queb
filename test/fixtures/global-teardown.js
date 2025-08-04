// test/fixtures/global-teardown.js - Global Jest teardown for clean exits

module.exports = async () => {
  // Force cleanup of any remaining resources
  const { execSync } = require('child_process');
  
  try {
    // Kill any test servers on common ports
    execSync('lsof -ti:8080 | xargs kill -9 2>/dev/null || true', { stdio: 'ignore' });
    execSync('lsof -ti:3000 | xargs kill -9 2>/dev/null || true', { stdio: 'ignore' });
  } catch (e) {
    // Ignore errors - ports may not be in use
  }
  
  // Clear any global timers
  if (global.testCleanupRegistry) {
    await global.testCleanupRegistry.cleanup();
  }
  
  // Force garbage collection if available
  if (global.gc) {
    global.gc();
  }
  
  console.log('✅ Global teardown complete');
};