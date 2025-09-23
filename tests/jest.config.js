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
      testMatch: [
        '**/suites/unit/unit.test.js',
        '**/suites/unit/check-env.test.js',
        '**/suites/unit/health-check.test.js',
        '**/suites/unit/backup.test.js',
        '**/suites/unit/restore.test.js',
        '**/suites/unit/setup-database.test.js'
      ],
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
    },
    {
      displayName: 'e2e',
      testMatch: ['**/suites/e2e/*.test.js'],
      testEnvironment: 'node',
      globals: {
        TextEncoder: TextEncoder,
        TextDecoder: TextDecoder
      },
      detectOpenHandles: true
    }
  ],
  forceExit: true,
  detectOpenHandles: true,
  testTimeout: 120000, // 2 minutes for e2e tests
  verbose: false
};


