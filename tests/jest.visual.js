module.exports = {
  testEnvironment: 'jsdom',
  testMatch: [
    '**/visual/**/*.js',
    '**/?(*.)+(visual|spec|test).js'
  ],
  setupFilesAfterEnv: ['<rootDir>/tests/visual-setup.js'],
  transform: {
    '^.+\\.js$': 'babel-jest',
    '^.+\\.jsx$': 'babel-jest'
  },
  moduleFileExtensions: ['js', 'jsx', 'json'],
  testTimeout: 30000,
  collectCoverage: false,
  // Visual regression testing setup
  setupFiles: ['jest-canvas-mock']
};