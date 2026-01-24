module.exports = {
  testEnvironment: 'node',
  testMatch: [
<<<<<<< HEAD
    '<rootDir>/tests/suites/unit/**/*.test.js',
    '<rootDir>/tests/unit/**/*.spec.js',
    '<rootDir>/tests/unit/**/*.test.js'
=======
    '**/__tests__/**/*.js',
    '**/?(*.)+(spec|test).js',
    '!**/integration/**',
    '!**/e2e/**',
    '!**/visual/**'
>>>>>>> 532fb1f (merged in website)
  ],
  collectCoverageFrom: [
    'src/**/*.js',
    '!src/**/*.test.js',
    '!src/**/*.spec.js'
  ],
  setupFilesAfterEnv: ['<rootDir>/tests/setup.js'],
  transform: {
    '^.+\\.js$': 'babel-jest'
  },
  moduleFileExtensions: ['js', 'json'],
  testTimeout: 10000
};