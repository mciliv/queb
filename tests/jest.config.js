module.exports = {
  rootDir: __dirname,
  testEnvironment: 'node',
  reporters: [
    'default',
    '<rootDir>/reporters/one-line-reporter.js'
  ],
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
  // Reduce noisy console during tests; reporter provides one-liners
  silent: true,
  transform: {
    '^.+\\.(js|jsx)$': ['babel-jest', { 
      presets: [
        ['@babel/preset-env', { 
          targets: { node: 'current' },
          modules: 'commonjs'
        }],
        ['@babel/preset-react', { runtime: 'automatic' }]
      ],
      plugins: ['@babel/plugin-transform-modules-commonjs']
    }]
  },
  // Ensure source files are transformed
  transformIgnorePatterns: [
    'node_modules/(?!(.*\\.(js|mjs)$))'
  ],
  // Include source files in transformation
  moduleNameMapper: {
    '^queb$': '<rootDir>/../package.json',
    '^@/(.*)$': '<rootDir>/../src/$1'
  },
  moduleFileExtensions: ['js', 'jsx', 'json'],
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
  // Fail on skipped tests and ensure all tests pass
  passWithNoTests: false,
  errorOnDeprecated: true,
  projects: [
    {
      displayName: 'unit-frontend',
      testMatch: ['**/suites/unit/camera*.test.js', '**/suites/unit/manual.test.js', '**/suites/unit/front-end*.test.js', '**/suites/unit/keyboard-shortcuts.test.js'],
      testEnvironment: 'jsdom',
      setupFilesAfterEnv: ['<rootDir>/support/fixtures/setup.js'],
      transform: {
        '^.+\\.(js|jsx)$': ['babel-jest', { 
          presets: [
            ['@babel/preset-env', { 
              targets: { node: 'current' },
              modules: 'commonjs'
            }],
            ['@babel/preset-react', { runtime: 'automatic' }]
          ],
          plugins: ['@babel/plugin-transform-modules-commonjs']
        }]
      },
      transformIgnorePatterns: [
        'node_modules/(?!(.*\\.(js|mjs)$))'
      ],
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
        '**/suites/unit/setup-database.test.js',
        '**/suites/unit/sdf-retriever.test.js'
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


