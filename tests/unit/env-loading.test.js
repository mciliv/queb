const fs = require('fs');
const path = require('path');

// Mock modules before requiring the services module
jest.mock('fs');
jest.mock('path', () => {
  const real = jest.requireActual('path');
  return {
    dirname: jest.fn((...a) => real.dirname(...a)),
    parse: jest.fn((...a) => real.parse(...a)),
    resolve: jest.fn((...a) => real.resolve(...a)),
    join: jest.fn((...a) => real.join(...a)),
    sep: real.sep,
    posix: real.posix,
  };
});

// Import the functions we want to test after mocking
// Note: We need to require the entire module to test the functions
const servicesModule = require('../../src/core/services');

// Extract the functions we want to test
const {
  findEnvFilePath,
  findProjectRoot,
  validateCriticalEnvironmentVariables,
  loadEnvironmentVariables
} = (() => {
  const mod = require('../../src/core/services');
  return {
    findEnvFilePath: mod.findEnvFilePath || (() => {}),
    findProjectRoot: mod.findProjectRoot || (() => {}),
    validateCriticalEnvironmentVariables: mod.validateCriticalEnvironmentVariables || (() => {}),
    loadEnvironmentVariables: mod.loadEnvironmentVariables || (() => {})
  };
})();

describe('Environment Variable Loading System', () => {
  let originalEnv;
  let consoleLogSpy;
  let consoleWarnSpy;
  let consoleErrorSpy;

  beforeEach(() => {
    // Save original environment
    originalEnv = { ...process.env };

    // Mock console methods
    consoleLogSpy = jest.spyOn(console, 'log').mockImplementation();
    consoleWarnSpy = jest.spyOn(console, 'warn').mockImplementation();
    consoleErrorSpy = jest.spyOn(console, 'error').mockImplementation();

    // Reset global flag
    delete global.__QUEB_ENV_LOADED__;
  });

  afterEach(() => {
    // Restore original environment
    process.env = originalEnv;

    // Restore console methods
    consoleLogSpy.mockRestore();
    consoleWarnSpy.mockRestore();
    consoleErrorSpy.mockRestore();

    // Clear mocks
    jest.clearAllMocks();
  });

  describe('findProjectRoot', () => {
    beforeEach(() => {
      // Reset path mocks to real implementations before each test
      const real = jest.requireActual('path');
      path.dirname.mockImplementation((p) => real.dirname(p));
      path.parse.mockReturnValue({ root: '/' });
      path.resolve.mockImplementation((...args) => real.resolve(...args));
    });

    it('should find project root by package.json', () => {
      fs.existsSync
        .mockReturnValueOnce(false) // first dir: no package.json
        .mockReturnValueOnce(false) // second dir: no package.json
        .mockReturnValueOnce(true);  // third dir: package.json found

      path.dirname
        .mockReturnValueOnce('/Users/project/src')
        .mockReturnValueOnce('/Users/project')
        .mockReturnValueOnce('/');

      const result = findProjectRoot();
      expect(result).toBe('/Users/project');
    });

    it('should fallback to hardcoded path if package.json not found', () => {
      fs.existsSync.mockReturnValue(false);
      const real = jest.requireActual('path');
      path.resolve.mockReturnValue('/fallback/project/root');

      const result = findProjectRoot();
      expect(path.resolve).toHaveBeenCalledWith(expect.any(String), '../..');
    });
  });

  describe('validateCriticalEnvironmentVariables', () => {
    it('should warn when OPENAI_API_KEY is missing', () => {
      delete process.env.OPENAI_API_KEY;
      process.env.OPENAI_MODEL = 'gpt-4';

      validateCriticalEnvironmentVariables();

      expect(consoleWarnSpy).toHaveBeenCalledWith(
        expect.stringContaining('Missing critical environment variables')
      );
      expect(consoleWarnSpy).toHaveBeenCalledWith(
        expect.stringContaining('OPENAI_API_KEY')
      );
    });

    it('should warn when OPENAI_MODEL is missing', () => {
      process.env.OPENAI_API_KEY = 'sk-test';
      delete process.env.OPENAI_MODEL;

      validateCriticalEnvironmentVariables();

      expect(consoleWarnSpy).toHaveBeenCalledWith(
        expect.stringContaining('Missing critical environment variables')
      );
      expect(consoleWarnSpy).toHaveBeenCalledWith(
        expect.stringContaining('OPENAI_MODEL')
      );
    });

    it('should log success when all critical vars are present', () => {
      process.env.OPENAI_API_KEY = 'sk-test';
      process.env.OPENAI_MODEL = 'gpt-4';

      validateCriticalEnvironmentVariables();

      expect(consoleLogSpy).toHaveBeenCalledWith(
        '✅ All critical environment variables loaded successfully'
      );
    });
  });
});
