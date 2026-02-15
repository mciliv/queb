/**
 * Global setup for all tests
 * This runs once before all test suites
 */

const fs = require('fs');
const path = require('path');

module.exports = async () => {
  // Ensure test directories exist
  const testDirs = [
    path.join(__dirname, '../coverage'),
    path.join(__dirname, '../reports'),
    path.join(__dirname, '../screenshots'),
    path.join(__dirname, '../../public/sdf_files/test')
  ];

  testDirs.forEach(dir => {
    if (!fs.existsSync(dir)) {
      fs.mkdirSync(dir, { recursive: true });
    }
  });

  // Set test environment variables
  process.env.NODE_ENV = 'test';
  process.env.LOG_LEVEL = process.env.LOG_LEVEL || 'error';
  process.env.TEST_MODE = 'true';

  // Initialize test database if needed
  if (process.env.TEST_DATABASE_URL) {
    console.log('Initializing test database...');
    // Database initialization logic here
  }

  console.log('Global test setup complete');
};
