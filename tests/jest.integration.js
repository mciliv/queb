const baseConfig = require('./jest.base');

module.exports = {
  ...baseConfig,
  
  // Display name for this project
  displayName: 'integration',
  
  // Test match patterns for integration tests
  testMatch: [
    '<rootDir>/tests/suites/integration/**/*.test.js',
    '<rootDir>/tests/integration/**/*.test.js'
  ],
  
  // Integration tests need more time
  testTimeout: 60000, // 60 seconds
  
  // Setup files specific to integration tests
  setupFilesAfterEnv: [
    ...baseConfig.setupFilesAfterEnv,
    '<rootDir>/tests/utils/integration-setup.js'
  ],
  
  // Coverage directory
  coverageDirectory: '<rootDir>/tests/coverage/integration',
  
  // Detect open handles for integration tests
  detectOpenHandles: true,
  
  // Integration tests might need to exit forcefully
  forceExit: true,
  
  // Run tests serially for integration tests to avoid conflicts
  maxWorkers: 1,
  
  // Integration test reporters
  reporters: (() => {
    const reporters = ['default'];
    // jest-junit reporter (optional - only if package is installed)
    try {
      require.resolve('jest-junit');
      reporters.push(['jest-junit', {
        outputDirectory: '<rootDir>/tests/reports/integration',
        outputName: 'junit.xml',
        classNameTemplate: '{classname}',
        titleTemplate: '{title}',
        ancestorSeparator: ' › ',
        usePathForSuiteName: true
      }]);
    } catch (e) {
      // jest-junit not installed, skip it
    }
    return reporters;
  })()
};
