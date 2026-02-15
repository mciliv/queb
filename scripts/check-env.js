#!/usr/bin/env node
/**
 * Check environment configuration and dependencies
 */

const net = require('net');

// Colored logging utilities
const log = {
  info: (msg) => console.log(`\x1b[32m✅ ${msg}\x1b[0m`),
  warn: (msg) => console.log(`\x1b[33m⚠️  ${msg}\x1b[0m`),
  error: (msg) => console.log(`\x1b[31m❌ ${msg}\x1b[0m`)
};

/**
 * Check if a port is in use
 * @param {number} port - Port to check
 * @returns {Promise<boolean>} - true if port is in use, false otherwise
 */
async function checkPort(port) {
  return new Promise((resolve) => {
    const server = net.createServer();
    
    server.once('error', (err) => {
      if (err.code === 'EADDRINUSE') {
        resolve(true); // Port is in use
      } else {
        resolve(false);
      }
    });
    
    server.once('listening', () => {
      server.close();
      resolve(false); // Port is available
    });
    
    server.listen(port, '127.0.0.1');
  });
}

/**
 * Check if a URL is accessible
 * @param {string} url - URL to check
 * @param {number} timeout - Timeout in milliseconds
 * @returns {Promise<boolean>} - true if URL is accessible, false otherwise
 */
async function checkUrl(url, timeout = 5000) {
  const controller = new AbortController();
  let timeoutId = null;
  try {
    timeoutId = setTimeout(() => controller.abort(), timeout);

    const response = await fetch(url, {
      signal: controller.signal,
      method: 'HEAD'
    });

    return response.ok;
  } catch (error) {
    return false;
  } finally {
    if (timeoutId) clearTimeout(timeoutId);
  }
}

// Main environment check
async function main() {
  log.info('Checking environment configuration...');
  
  const port = parseInt(process.env.PORT) || 8080;
  const isPortInUse = await checkPort(port);
  
  if (isPortInUse) {
    log.warn(`Port ${port} is already in use`);
  } else {
    log.info(`Port ${port} is available`);
  }
  
  // Check if we're in production or development
  const env = process.env.NODE_ENV || 'development';
  log.info(`Environment: ${env}`);
  
  log.info('Environment check complete');
}

if (require.main === module) {
  main().catch(error => {
    log.error(`Environment check failed: ${error.message}`);
    process.exit(1);
  });
}

module.exports = { log, checkPort, checkUrl };

