module.exports = {
  rootDir: __dirname,
  testEnvironment: 'jsdom',
  setupFilesAfterEnv: ['<rootDir>/support/fixtures/setup.js'],
  globalTeardown: '<rootDir>/support/fixtures/global-teardown.js',
  testMatch: [
    '<rootDir>/**/*.test.js'
  ],
  testPathIgnorePatterns: [
    '<rootDir>/../dist/',
    '<rootDir>/../node_modules/'
  ],
  collectCoverageFrom: [
    '<rootDir>/../backend/**/*.js',
    '<rootDir>/../frontend/**/*.js',
    '!<rootDir>/../**/node_modules/**',
    '!<rootDir>/../**/dist/**'
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
    '^mol$': '<rootDir>/../package.json'
  },
  projects: [
    {
      displayName: 'unit-frontend',
      testMatch: ['**/suites/unit/camera*.test.js', '**/suites/unit/manual.test.js', '**/suites/unit/front-end*.test.js', '**/suites/unit/keyboard-shortcuts.test.js'],
      testEnvironment: 'jsdom',
      setupFilesAfterEnv: ['<rootDir>/support/fixtures/setup.js'],
      globals: {
        TextEncoder: TextEncoder,
        TextDecoder: TextDecoder
      },
      detectOpenHandles: true
    },
    {
      displayName: 'unit-backend',
      testMatch: ['**/suites/unit/unit.test.js'],
      testEnvironment: 'node',
      globals: {
        TextEncoder: TextEncoder,
        TextDecoder: TextDecoder
      },
      detectOpenHandles: true
    },
    {
      displayName: 'integration',
      testMatch: ['**/suites/integration/*.test.js'],
      testEnvironment: 'node',
      globals: {
        TextEncoder: TextEncoder,
        TextDecoder: TextDecoder
      },
      detectOpenHandles: true
    },
    {
      displayName: 'smoke',
      testMatch: ['**/smoke.test.js'],
      testEnvironment: 'node',
      detectOpenHandles: true
    }
  ],
  forceExit: true,
  detectOpenHandles: true,
  testTimeout: 30000,
  verbose: false
};


