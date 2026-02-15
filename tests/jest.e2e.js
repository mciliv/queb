module.exports = {
  testEnvironment: 'node',
  testMatch: [
    '**/e2e/**/*.js',
    '**/?(*.)+(e2e|spec|test).js'
  ],
  setupFilesAfterEnv: ['<rootDir>/tests/setup.js'],
  transform: {
    '^.+\\.js$': 'babel-jest'
  },
  moduleFileExtensions: ['js', 'json'],
  testTimeout: 60000, // E2E tests can take longer
  collectCoverage: false,
  globalSetup: '<rootDir>/tests/e2e-setup.js',
  globalTeardown: '<rootDir>/tests/e2e-teardown.js'
};