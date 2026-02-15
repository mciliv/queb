/**
 * Integration test setup
 * Additional setup for integration tests
 */

const path = require('path');
const { spawn } = require('child_process');

// Store server process
let serverProcess = null;

beforeAll(async () => {
  // Start test server if not already running
  const serverPort = process.env.TEST_SERVER_PORT || 8081;
  
  if (!process.env.USE_EXISTING_SERVER) {
    console.log(`Starting test server on port ${serverPort}...`);
    
    serverProcess = spawn('node', [
      path.join(__dirname, '../../src/server/api/server.js')
    ], {
      env: {
        ...process.env,
        PORT: serverPort,
        NODE_ENV: 'test',
        LOG_LEVEL: 'error'
      },
      stdio: 'pipe'
    });

    // Wait for server to start
    await new Promise((resolve, reject) => {
      let output = '';
      
      const checkServer = async () => {
        try {
          const response = await fetch(`http://localhost:${serverPort}/health`);
          if (response.ok) {
            resolve();
          } else {
            setTimeout(checkServer, 500);
          }
        } catch (error) {
          setTimeout(checkServer, 500);
        }
      };
      
      serverProcess.stdout.on('data', (data) => {
        output += data.toString();
        if (output.includes('Server started') || output.includes('listening')) {
          setTimeout(checkServer, 1000);
        }
      });
      
      serverProcess.stderr.on('data', (data) => {
        console.error('Server error:', data.toString());
      });
      
      serverProcess.on('error', reject);
      
      // Start checking after a delay
      setTimeout(checkServer, 2000);
      
      // Timeout after 30 seconds
      setTimeout(() => reject(new Error('Server failed to start')), 30000);
    });
    
    console.log('Test server started successfully');
  }
  
  // Set base URL for tests
  global.TEST_BASE_URL = `http://localhost:${serverPort}`;
});

afterAll(async () => {
  // Stop test server if we started it
  if (serverProcess) {
    console.log('Stopping test server...');
    
    serverProcess.kill('SIGTERM');
    
    // Wait for process to exit
    await new Promise((resolve) => {
      serverProcess.on('exit', resolve);
      setTimeout(resolve, 5000); // Force resolve after 5 seconds
    });
    
    console.log('Test server stopped');
  }
});

// Helper to create authenticated requests
global.createAuthenticatedRequest = (endpoint, options = {}) => {
  return fetch(`${global.TEST_BASE_URL}${endpoint}`, {
    ...options,
    headers: {
      'Content-Type': 'application/json',
      'Authorization': `Bearer ${process.env.TEST_API_KEY || 'test-key'}`,
      ...options.headers
    }
  });
};

// Helper to wait for async operations
global.waitForResponse = async (promiseFn, expectedStatus = 200, retries = 3) => {
  for (let i = 0; i < retries; i++) {
    try {
      const response = await promiseFn();
      if (response.status === expectedStatus) {
        return response;
      }
    } catch (error) {
      if (i === retries - 1) throw error;
      await new Promise(resolve => setTimeout(resolve, 1000));
    }
  }
  throw new Error(`Failed to get expected status ${expectedStatus} after ${retries} attempts`);
};
