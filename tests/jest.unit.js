const baseConfig = require('./jest.base');

module.exports = {
  ...baseConfig,
  
  // Display name for this project
  displayName: 'unit',
  
  // Test match patterns for unit tests
  testMatch: [
    '<rootDir>/tests/suites/unit/**/*.test.js',
    '<rootDir>/tests/unit/**/*.spec.js'
  ],
  
  // Unit tests should be fast
  testTimeout: 10000, // 10 seconds
  
  // Coverage disabled for faster test runs (enable when needed)
  collectCoverage: false,
  
  // Coverage directory
  coverageDirectory: '<rootDir>/tests/coverage/unit',
  
  // Don't detect open handles for unit tests (they should be fast)
  detectOpenHandles: false,
  
  // Unit tests specific reporters
  reporters: (() => {
    const reporters = ['default'];
    // jest-junit reporter (optional - only if package is installed)
    try {
      require.resolve('jest-junit');
      reporters.push(['jest-junit', {
        outputDirectory: '<rootDir>/tests/reports/unit',
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
