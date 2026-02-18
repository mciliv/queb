/**
 * Global setup for E2E tests
 */

module.exports = async () => {
  // Start test database or services
  console.log('Setting up E2E test environment...');

  // Set test-specific environment variables
  process.env.NODE_ENV = 'test';
  process.env.E2E_TEST = 'true';

  // Could start a test server here if needed
  // global.testServer = await startTestServer();
};