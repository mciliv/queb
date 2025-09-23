#!/usr/bin/env node

/**
 * Quick Environment Check - Node.js version
 * Replaces scripts/check-env.sh
 */

const { execSync } = require('child_process');
const net = require('net');

// Colored logging utilities (replacing bash utils.sh)
const log = {
  info: (msg) => console.log(`\x1b[32m✅ ${msg}\x1b[0m`),
  warn: (msg) => console.log(`\x1b[33m⚠️  ${msg}\x1b[0m`),
  error: (msg) => console.log(`\x1b[31m❌ ${msg}\x1b[0m`)
};

async function checkPort(port) {
  return new Promise((resolve) => {
    const server = net.createServer();
    server.listen(port, '127.0.0.1', () => {
      server.close();
      resolve(false); // Port is free
    });
    server.on('error', () => resolve(true)); // Port is in use
  });
}

async function checkUrl(url, timeout = 5000) {
  try {
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), timeout);
    
    const response = await fetch(url, { 
      signal: controller.signal,
      method: 'HEAD'
    });
    
    clearTimeout(timeoutId);
    return response.ok;
  } catch (error) {
    return false;
  }
}

async function main() {
  console.log('🔍 Environment Status');
  console.log('====================');

  // Environment
  const nodeEnv = process.env.NODE_ENV;
  if (nodeEnv === 'production') {
    console.log('🚨 PRODUCTION');
  } else {
    console.log('💻 DEVELOPMENT');
  }

  // Local server
  const serverRunning = await checkPort(8080);
  if (serverRunning) {
    console.log('🟢 Server: RUNNING (http://localhost:8080)');
  } else {
    console.log('🔴 Server: NOT RUNNING (run \'d\' to start)');
  }

  // Production check
  const prodAccessible = await checkUrl('https://queb.space');
  if (prodAccessible) {
    console.log('🟢 Production: ACCESSIBLE (https://queb.space)');
  } else {
    console.log('🟡 Production: NOT ACCESSIBLE');
  }

  // Git status
  try {
    if (require('fs').existsSync('.git')) {
      const changes = execSync('git status --porcelain', { encoding: 'utf8' })
        .split('\n')
        .filter(line => line.trim()).length;
      
      if (changes === 0) {
        console.log('✅ Git: Clean');
      } else {
        console.log(`📝 Git: ${changes} changes`);
      }
    }
  } catch (error) {
    // Ignore git errors
  }

  console.log('');
  console.log('Quick commands: d (dev) | d deploy (deploy) | d clean (cleanup)');
}

if (require.main === module) {
  main().catch(console.error);
}

module.exports = { log, checkPort, checkUrl };
