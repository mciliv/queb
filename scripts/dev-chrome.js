#!/usr/bin/env node

/**
 * Chrome Auto-Reload Development Script
 * 
 * This script starts the development server with Chrome auto-reload integration.
 * It replaces the live-reload system with Chrome DevTools Protocol integration.
 */

const { spawn } = require('child_process');
const path = require('path');
const DevServer = require('../src/core/DevServer');
const configuration = require('../src/core/Configuration');

class ChromeDevEnvironment {
  constructor() {
    this.serverProcess = null;
    this.devServer = null;
    this.isShuttingDown = false;
  }

  async start() {
    console.log('🚀 Starting Chrome Auto-Reload Development Environment...');
    
    try {
      // Start the application server
      await this.startAppServer();
      
      // Wait for server to be ready
      await this.waitForServer();
      
      // Start Chrome auto-reload
      await this.startChromeAutoReload();
      
      // Setup shutdown handlers
      this.setupShutdownHandlers();
      
      console.log('\n✅ Development environment ready!');
      console.log('📱 Application: http://localhost:8080');
      console.log('🔒 HTTPS: https://localhost:3001');
      console.log('🌐 Chrome will auto-reload on file changes');
      console.log('\nPress Ctrl+C to stop\n');
      
    } catch (error) {
      console.error('❌ Failed to start development environment:', error.message);
      await this.shutdown();
      process.exit(1);
    }
  }

  async startAppServer() {
    console.log('🔧 Starting application server...');
    
    const serverPath = path.resolve(__dirname, '../src/server/api/server.js');
    
    this.serverProcess = spawn('node', [serverPath], {
      stdio: ['inherit', 'inherit', 'inherit'],
      env: {
        ...process.env,
        NODE_ENV: 'development'
      }
    });
    
    this.serverProcess.on('error', (error) => {
      console.error('❌ Server process error:', error);
    });
    
    this.serverProcess.on('exit', (code, signal) => {
      if (!this.isShuttingDown) {
        console.error(`❌ Server exited unexpectedly (code: ${code}, signal: ${signal})`);
        this.shutdown();
      }
    });
    
    console.log('✅ Application server started');
  }

  async waitForServer(maxAttempts = 30) {
    console.log('⏳ Waiting for server to be ready...');
    
    for (let i = 0; i < maxAttempts; i++) {
      try {
        const response = await fetch('http://localhost:8080/health');
        if (response.ok) {
          console.log('✅ Server is ready');
          return;
        }
      } catch (error) {
        // Server not ready yet
      }
      
      await new Promise(resolve => setTimeout(resolve, 1000));
    }
    
    throw new Error('Server failed to start within timeout period');
  }

  async startChromeAutoReload() {
    console.log('🌐 Starting Chrome auto-reload...');
    
    this.devServer = new DevServer({
      appUrl: 'http://localhost:8080',
      chromePort: 9222,
      debounceMs: 1000, // Longer debounce to avoid rapid reloads
      watchDirs: [
        'src/client',
        'src/core',
        'src/server'
      ],
      ignorePatterns: [
        /node_modules/,
        /\.git/,
        /dist/,
        /logs/,
        /temp/,
        /\.DS_Store/,
        /\.map$/,
        /package-lock\.json$/
      ]
    });
    
    // Handle DevServer events
    this.devServer.on('ready', () => {
      console.log('✅ Chrome auto-reload ready');
    });
    
    this.devServer.on('error', (error) => {
      console.warn('⚠️ Chrome auto-reload warning:', error.message);
    });
    
    await this.devServer.start();
  }

  setupShutdownHandlers() {
    const signals = ['SIGINT', 'SIGTERM', 'SIGUSR2'];
    
    signals.forEach(signal => {
      process.on(signal, async () => {
        console.log(`\n🛑 Received ${signal}, shutting down gracefully...`);
        await this.shutdown();
        process.exit(0);
      });
    });
  }

  async shutdown() {
    if (this.isShuttingDown) return;
    this.isShuttingDown = true;
    
    console.log('🧹 Cleaning up development environment...');
    
    // Stop Chrome auto-reload
    if (this.devServer) {
      try {
        await this.devServer.stop();
        console.log('✅ Chrome auto-reload stopped');
      } catch (error) {
        console.warn('⚠️ Error stopping DevServer:', error.message);
      }
    }
    
    // Stop application server
    if (this.serverProcess) {
      try {
        this.serverProcess.kill('SIGTERM');
        
        // Give it time to shutdown gracefully
        await new Promise(resolve => {
          const timeout = setTimeout(() => {
            this.serverProcess.kill('SIGKILL');
            resolve();
          }, 5000);
          
          this.serverProcess.on('exit', () => {
            clearTimeout(timeout);
            resolve();
          });
        });
        
        console.log('✅ Application server stopped');
      } catch (error) {
        console.warn('⚠️ Error stopping server:', error.message);
      }
    }
    
    console.log('✅ Development environment stopped');
  }
}

// Start the development environment
if (require.main === module) {
  const env = new ChromeDevEnvironment();
  env.start().catch(error => {
    console.error('💥 Fatal error:', error);
    process.exit(1);
  });
}

module.exports = ChromeDevEnvironment;
