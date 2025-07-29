module.exports = {
  testEnvironment: 'jsdom',
  setupFilesAfterEnv: ['<rootDir>/test/fixtures/setup.js'],
  testMatch: [
    '**/test/**/*.test.js',
    '!**/node_modules/**'
  ],
  collectCoverageFrom: [
    'backend/**/*.js',
    'frontend/**/*.js',
    '!**/node_modules/**',
    '!**/coverage/**'
  ],
  testTimeout: 30000,
  verbose: true,
  detectOpenHandles: true,
  testEnvironmentOptions: {
    url: 'http://localhost:8080'
  },
  transform: {
    '^.+\\.(js|jsx)$': ['babel-jest']
  },
  transformIgnorePatterns: [
    'node_modules/(?!(.*\\.(js|mjs)$))'
  ],
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
      testMatch: ['**/test/unit/camera*.test.js', '**/test/unit/manual.test.js', '**/test/unit/front-end*.test.js'],
      testEnvironment: 'jsdom',
      setupFilesAfterEnv: ['<rootDir>/test/fixtures/setup.js'],
          detectOpenHandles: true
    },
    {
      displayName: 'unit-backend', 
      testMatch: ['**/test/unit/unit.test.js'],
      testEnvironment: 'node',
      globals: {
        TextEncoder: TextEncoder,
        TextDecoder: TextDecoder
      },
          detectOpenHandles: true
    },
    {
      displayName: 'integration',
      testMatch: ['**/test/integration/*.test.js'],
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
  ]
}; 