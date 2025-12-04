module.exports = {
  // Root directory for tests
  rootDir: '..',
  
  // Test environment defaults
  testEnvironment: 'node',
  
  // Transform files with babel-jest
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
  
  // Module resolution
  moduleNameMapper: {
    '^@/(.*)$': '<rootDir>/src/$1',
    '^@test/(.*)$': '<rootDir>/tests/utils/$1',
    '^@fixtures/(.*)$': '<rootDir>/tests/fixtures/$1'
  },
  
  // File extensions to consider
  moduleFileExtensions: ['js', 'jsx', 'json'],
  
  // Ignore patterns
  testPathIgnorePatterns: [
    '/node_modules/',
    '/dist/',
    '/build/',
    '/.cache/'
  ],
  
  // Transform ignore patterns
  transformIgnorePatterns: [
    'node_modules/(?!(@babel|@testing-library)/)'
  ],
  
  // Coverage collection
  collectCoverageFrom: [
    'src/**/*.{js,jsx}',
    '!src/**/*.test.{js,jsx}',
    '!src/**/*.spec.{js,jsx}',
    '!src/**/test/**',
    '!src/**/tests/**',
    '!src/**/__tests__/**',
    '!src/**/node_modules/**'
  ],
  
  // Coverage thresholds (disabled - enable when coverage improves)
  // coverageThreshold: {
  //   global: {
  //     branches: 70,
  //     functions: 70,
  //     lines: 70,
  //     statements: 70
  //   }
  // },
  
  // Global setup/teardown
  globalSetup: '<rootDir>/tests/utils/setup.js',
  globalTeardown: '<rootDir>/tests/utils/teardown.js',
  
  // Setup files after env
  setupFilesAfterEnv: ['<rootDir>/tests/utils/setup-test-env.js'],
  
  // Globals
  globals: {
    TextEncoder: TextEncoder,
    TextDecoder: TextDecoder
  },
  
  // Error handling
  errorOnDeprecated: true,
  
  // Timeouts
  testTimeout: 30000, // 30 seconds default
  
  // Verbose output
  verbose: true,
  
  // Detect open handles (useful for debugging)
  detectOpenHandles: false,
  
  // Force exit after tests complete
  forceExit: false,
  
  // Clear mocks automatically
  clearMocks: true,
  
  // Restore mocks automatically
  restoreMocks: true,
  
  // Reset mocks automatically
  resetMocks: true
};
