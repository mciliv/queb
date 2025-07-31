module.exports = {
  testEnvironment: 'jsdom',
  setupFilesAfterEnv: ['<rootDir>/test/fixtures/setup.js'],
  globalTeardown: '<rootDir>/test/fixtures/global-teardown.js',
  testMatch: [
    '<rootDir>/test/**/*.test.js'
  ],
  moduleDirectories: ['node_modules'],
  testPathIgnorePatterns: ['/node_modules/'],
  collectCoverageFrom: [
    'backend/**/*.js',
    'frontend/**/*.js',
    '!**/node_modules/**'
  ],
  transform: {
    '^.+\\.(js|jsx)$': ['babel-jest']
  },
  transformIgnorePatterns: [
    'node_modules/(?!(.*\\.(js|mjs)$))'
  ],
  testEnvironmentOptions: {
    url: 'http://localhost:8080'
  },
  globals: {
    TextEncoder: TextEncoder,
    TextDecoder: TextDecoder
  },
  forceExit: true,
  detectOpenHandles: true,
  testTimeout: 30000
}; 