// test/fixtures/global-teardown.js - Global Jest teardown for clean exits

module.exports = async () => {
  // Force cleanup of any remaining resources
  const { execSync } = require('child_process');
  const fs = require('fs');
  const path = require('path');
  const ROOT = path.resolve(__dirname, '../../..');
  const TMP_DIR = path.join(ROOT, '.tmp');

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

  // Remove any Chrome userDataDir left by tests
  try {
    if (fs.existsSync(TMP_DIR)) {
      for (const name of fs.readdirSync(TMP_DIR)) {
        if (
          name.startsWith('chrome-molecular-profile-') ||
          name === 'seamless-chrome-profile' ||
          name === 'single-tab-chrome'
        ) {
          const fullPath = path.join(TMP_DIR, name);
          try { fs.rmSync(fullPath, { recursive: true, force: true }); } catch (_) {}
        }
      }
    }
  } catch (_) {}

  // Force garbage collection if available
  if (global.gc) {
    global.gc();
  }
  
  console.log('✅ Global teardown complete');
};