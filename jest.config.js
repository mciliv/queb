module.exports = {
  testEnvironment: 'jsdom',
  setupFilesAfterEnv: ['<rootDir>/test/support/fixtures/setup.js'],
  globalTeardown: '<rootDir>/test/support/fixtures/global-teardown.js',
  testMatch: [
    '<rootDir>/test/**/*.test.js'
  ],
  testPathIgnorePatterns: [
    '<rootDir>/dist/',
    '<rootDir>/node_modules/'
  ],
  collectCoverageFrom: [
    'backend/**/*.js',
    'frontend/**/*.js',
    '!**/node_modules/**',
    '!**/dist/**'
  ],
  transform: {
    '^.+\\.(js|jsx)$': ['babel-jest']
  },
  transformIgnorePatterns: [
    'node_modules/(?!(.*\\.(js|mjs)$))',
    '!frontend/'
  ],
  modulePathIgnorePatterns: [
    '<rootDir>/dist/'
  ],
  testEnvironmentOptions: {
    url: 'http://localhost:8080'
  },
  globals: {
    TextEncoder: TextEncoder,
    TextDecoder: TextDecoder
  },
  moduleNameMapper: {
    '^mol$': '<rootDir>/package.json'
  },
  projects: [
    {
      displayName: 'unit-frontend',
      testMatch: ['**/test/suites/unit/camera*.test.js', '**/test/suites/unit/manual.test.js', '**/test/suites/unit/front-end*.test.js'],
      testEnvironment: 'jsdom',
      setupFilesAfterEnv: ['<rootDir>/test/support/fixtures/setup.js'],
      globals: {
        TextEncoder: TextEncoder,
        TextDecoder: TextDecoder
      },
      detectOpenHandles: true
    },
    {
      displayName: 'unit-backend', 
      testMatch: ['**/test/suites/unit/unit.test.js'],
      testEnvironment: 'node',
      globals: {
        TextEncoder: TextEncoder,
        TextDecoder: TextDecoder
      },
      detectOpenHandles: true
    },
    {
      displayName: 'integration',
      testMatch: ['**/test/suites/integration/*.test.js'],
      testEnvironment: 'node',
      globals: {
        TextEncoder: TextEncoder,
        TextDecoder: TextDecoder
      },
      detectOpenHandles: true
    },
    {
      displayName: 'smoke',
      testMatch: ['**/test/smoke.test.js'],
      testEnvironment: 'node',
      detectOpenHandles: true
    }
  ],
  forceExit: true,
  detectOpenHandles: true,
  testTimeout: 30000,
  verbose: false
}; 