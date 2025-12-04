const baseConfig = require('./jest.base');

module.exports = {
  ...baseConfig,
  
  // Display name for this project
  displayName: 'visual',
  
  // Test match patterns for visual tests
  testMatch: [
    '<rootDir>/tests/suites/visual/**/*.test.js',
    '<rootDir>/tests/suites/visual/**/*.visual.spec.js'
  ],
  
  // Visual tests need time for rendering
  testTimeout: 90000, // 90 seconds
  
  // Setup files specific to visual tests
  setupFilesAfterEnv: [
    ...baseConfig.setupFilesAfterEnv,
    '<rootDir>/tests/utils/visual-setup.js'
  ],
  
  // Coverage directory
  coverageDirectory: '<rootDir>/tests/coverage/visual',
  
  // Don't collect coverage for visual tests
  collectCoverage: false,
  
  // Run visual tests serially to avoid screenshot conflicts
  maxWorkers: 1,
  
  // Visual test reporters
  reporters: [
    'default',
    ['jest-junit', {
      outputDirectory: '<rootDir>/tests/reports/visual',
      outputName: 'junit.xml',
      classNameTemplate: '{classname}',
      titleTemplate: '{title}',
      ancestorSeparator: ' › ',
      usePathForSuiteName: true
    }]
  ],
  
  // Global variables for visual tests
  globals: {
    ...baseConfig.globals,
    VISUAL_REGRESSION_THRESHOLD: process.env.VISUAL_REGRESSION_THRESHOLD || 0.01,
    UPDATE_SNAPSHOTS: process.env.UPDATE_SNAPSHOTS === 'true'
  }
};
