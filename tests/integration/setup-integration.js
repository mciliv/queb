/**
 * Integration Test Setup
 * Handles server startup and teardown for integration tests
 */

const { spawn } = require('child_process');
const path = require('path');

class IntegrationTestSetup {
  constructor() {
    this.serverProcess = null;
    this.baseUrl = 'http://localhost:8080';
  }

  async startServer() {
    return new Promise((resolve, reject) => {
      
      console.log('Starting server process...');
      this.serverProcess = spawn('node', ['backend/api/server.js'], {
        cwd: path.join(__dirname, '..', '..'),
        env: { 
          ...process.env, 
          NODE_ENV: 'test',
          INTEGRATION_TEST: 'true',
          PORT: '8080',
          NODE_OPTIONS: ''  // Disable debugger
        },
        stdio: ['pipe', 'pipe', 'pipe']
      });
      
      console.log('Server process started with PID:', this.serverProcess.pid);

      let serverReady = false;

      this.serverProcess.stdout.on('data', (data) => {
        const output = data.toString();
        console.log('[Server Output]', output.trim());
        
        if (output.includes('Server running on port') || output.includes('listening on')) {
          if (!serverReady) {
            serverReady = true;
            console.log('✅ Server ready for integration tests');
            resolve();
          }
        }
      });

      this.serverProcess.stderr.on('data', (data) => {
        const error = data.toString().trim();
        console.error(`[Server Stderr] ${error}`);
      });

      this.serverProcess.on('error', (error) => {
        console.error('❌ Failed to start server:', error);
        reject(error);
      });

      this.serverProcess.on('exit', (code) => {
        console.log(`Server process exited with code ${code}`);
        if (code !== 0 && !serverReady) {
          reject(new Error(`Server exited with code ${code}`));
        }
      });

      // Timeout after 30 seconds
      setTimeout(() => {
        if (!serverReady) {
          this.stopServer();
          reject(new Error('Server start timeout'));
        }
      }, 30000);
    });
  }

  stopServer() {
    if (this.serverProcess) {
      this.serverProcess.kill('SIGTERM');
      this.serverProcess = null;
    }
  }

  async waitForServer() {
    const maxAttempts = 30;
    let attempts = 0;
    
    while (attempts < maxAttempts) {
      try {
        const response = await fetch(`${this.baseUrl}/health`);
        if (response.ok) {
          return true;
        }
      } catch (error) {
        // Server not ready yet
      }
      
      await new Promise(resolve => setTimeout(resolve, 1000));
      attempts++;
    }
    
    throw new Error('Server did not become ready in time');
  }
}

module.exports = IntegrationTestSetup;