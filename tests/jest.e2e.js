const baseConfig = require('./jest.base');

module.exports = {
  ...baseConfig,
  
  // Display name for this project
  displayName: 'e2e',
  
  // Test match patterns for e2e tests
  testMatch: [
    '<rootDir>/tests/suites/e2e/**/*.test.js',
    '<rootDir>/tests/suites/e2e/**/*.e2e.spec.js'
  ],
  
  // E2E tests need the most time
  testTimeout: 120000, // 2 minutes
  
  // E2E tests should run in jsdom environment
  testEnvironment: 'node',
  
  // Setup files specific to e2e tests
  setupFilesAfterEnv: [
    ...baseConfig.setupFilesAfterEnv,
    '<rootDir>/tests/utils/e2e-setup.js'
  ],
  
  // Coverage directory
  coverageDirectory: '<rootDir>/tests/coverage/e2e',
  
  // Don't collect coverage for e2e tests by default
  collectCoverage: false,
  
  // Detect open handles
  detectOpenHandles: true,
  
  // E2E tests should exit forcefully
  forceExit: true,
  
  // Run e2e tests serially
  maxWorkers: 1,
  
  // E2E test reporters
  reporters: [
    'default',
    ['jest-junit', {
      outputDirectory: '<rootDir>/tests/reports/e2e',
      outputName: 'junit.xml',
      classNameTemplate: '{classname}',
      titleTemplate: '{title}',
      ancestorSeparator: ' › ',
      usePathForSuiteName: true
    }]
  ],
  
  // Global variables for e2e tests
  globals: {
    ...baseConfig.globals,
    BASE_URL: process.env.BASE_URL || 'http://localhost:8080',
    HEADLESS: process.env.HEADLESS !== 'false'
  }
};
