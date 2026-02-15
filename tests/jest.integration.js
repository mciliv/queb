module.exports = {
  testEnvironment: 'node',
  testMatch: [
    '**/integration/**/*.js',
    '**/?(*.)+(integration|spec|test).js'
  ],
  setupFilesAfterEnv: ['<rootDir>/tests/setup.js'],
  transform: {
    '^.+\\.js$': 'babel-jest'
  },
  moduleFileExtensions: ['js', 'json'],
  testTimeout: 30000, // Longer timeout for integration tests
  collectCoverage: false // Integration tests often have different coverage needs
};