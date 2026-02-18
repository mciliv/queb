const fs = require('fs');
const path = require('path');

// Mock modules before requiring the services module
jest.mock('fs');
jest.mock('dotenv');
jest.mock('path');

// Import the functions we want to test after mocking
// Note: We need to require the entire module to test the functions
const servicesModule = require('../../src/core/services');

// Extract the functions we want to test
const {
  detectCloudEnvironment,
  findEnvFilePath,
  findProjectRoot,
  validateCriticalEnvironmentVariables,
  loadEnvironmentVariables
} = (() => {
  // Re-require the module to get access to the functions
  // This is a bit of a hack but necessary for testing internal functions
  const mod = require('../../src/core/services');
  return {
    detectCloudEnvironment: mod.detectCloudEnvironment || (() => {}),
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

  describe('detectCloudEnvironment', () => {
    it('should detect Google App Engine', () => {
      process.env.GAE_APPLICATION = 'my-app';
      expect(detectCloudEnvironment()).toBe(true);
    });

    it('should detect Google Cloud Project', () => {
      process.env.GOOGLE_CLOUD_PROJECT = 'my-project';
      expect(detectCloudEnvironment()).toBe(true);
    });

    it('should detect Cloud Run', () => {
      process.env.K_SERVICE = 'my-service';
      expect(detectCloudEnvironment()).toBe(true);
    });

    it('should detect Cloud Functions', () => {
      process.env.FUNCTION_NAME = 'my-function';
      expect(detectCloudEnvironment()).toBe(true);
    });

    it('should detect AWS Lambda', () => {
      process.env.AWS_EXECUTION_ENV = 'AWS_Lambda_nodejs';
      expect(detectCloudEnvironment()).toBe(true);
    });

    it('should detect Vercel', () => {
      process.env.VERCEL = '1';
      expect(detectCloudEnvironment()).toBe(true);
    });

    it('should detect Netlify', () => {
      process.env.NETLIFY = 'true';
      expect(detectCloudEnvironment()).toBe(true);
    });

    it('should detect production NODE_ENV', () => {
      process.env.NODE_ENV = 'production';
      expect(detectCloudEnvironment()).toBe(true);
    });

    it('should not detect cloud environment for local development', () => {
      process.env.NODE_ENV = 'development';
      expect(detectCloudEnvironment()).toBe(false);
    });

    it('should not detect cloud environment when no cloud vars are set', () => {
      expect(detectCloudEnvironment()).toBe(false);
    });
  });

  describe('findProjectRoot', () => {
    beforeEach(() => {
      // Reset path mocks
      path.dirname.mockImplementation((p) => p.replace(/\/[^/]*$/, '') || '/');
      path.parse.mockReturnValue({ root: '/' });
      path.resolve.mockImplementation((...args) => args.join('/'));
    });

    it('should find project root by package.json', () => {
      const mockDir = '/Users/project/src/core';
      fs.existsSync
        .mockReturnValueOnce(false) // /Users/project/src/core/package.json
        .mockReturnValueOnce(false) // /Users/project/src/package.json
        .mockReturnValueOnce(true);  // /Users/project/package.json

      path.dirname
        .mockReturnValueOnce('/Users/project/src')     // from src/core
        .mockReturnValueOnce('/Users/project')         // from src
        .mockReturnValueOnce('/');                     // from project

      // Mock __dirname to be src/core
      Object.defineProperty(require('path'), '__dirname', {
        get: () => mockDir
      });

      const result = findProjectRoot();
      expect(result).toBe('/Users/project');
    });

    it('should fallback to hardcoded path if package.json not found', () => {
      fs.existsSync.mockReturnValue(false);
      path.resolve.mockReturnValue('/fallback/project/root');

      const result = findProjectRoot();
      expect(path.resolve).toHaveBeenCalledWith(expect.any(String), '../..');
    });
  });

  describe('findEnvFilePath', () => {
    it('should return .env path relative to project root', () => {
      const mockProjectRoot = '/Users/project';
      const mockEnvPath = '/Users/project/.env';

      // Mock findProjectRoot to return our test path
      const originalFindProjectRoot = servicesModule.findProjectRoot;
      servicesModule.findProjectRoot = jest.fn().mockReturnValue(mockProjectRoot);
      path.resolve.mockReturnValue(mockEnvPath);

      const result = findEnvFilePath();
      expect(result).toBe(mockEnvPath);
      expect(path.resolve).toHaveBeenCalledWith(mockProjectRoot, '.env');
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

  describe('loadEnvironmentVariables', () => {
    let dotenvMock;

    beforeEach(() => {
      dotenvMock = require('dotenv');
      dotenvMock.config.mockReturnValue({});

      // Mock the helper functions
      servicesModule.detectCloudEnvironment = jest.fn();
      servicesModule.findEnvFilePath = jest.fn();
      servicesModule.validateCriticalEnvironmentVariables = jest.fn();
    });

    it('should log cloud environment message and return early', () => {
      servicesModule.detectCloudEnvironment.mockReturnValue(true);

      loadEnvironmentVariables();

      expect(consoleLogSpy).toHaveBeenCalledWith(
        '🏭 Running in production/cloud environment - using SHELL ENVIRONMENT variables only'
      );
      expect(servicesModule.findEnvFilePath).not.toHaveBeenCalled();
      expect(dotenvMock.config).not.toHaveBeenCalled();
    });

    it('should load .env file for local development', () => {
      servicesModule.detectCloudEnvironment.mockReturnValue(false);
      servicesModule.findEnvFilePath.mockReturnValue('/path/to/.env');
      fs.existsSync.mockReturnValue(true);

      loadEnvironmentVariables();

      expect(consoleLogSpy).toHaveBeenCalledWith(
        '🏠 Running in local development - SHELL ENVIRONMENT takes priority, .env file provides fallbacks'
      );
      expect(consoleLogSpy).toHaveBeenCalledWith(
        '📄 Loading environment variables from: /path/to/.env'
      );
      expect(dotenvMock.config).toHaveBeenCalledWith({
        path: '/path/to/.env',
        override: false
      });
      expect(servicesModule.validateCriticalEnvironmentVariables).toHaveBeenCalled();
    });

    it('should warn when .env file does not exist', () => {
      servicesModule.detectCloudEnvironment.mockReturnValue(false);
      servicesModule.findEnvFilePath.mockReturnValue('/path/to/.env');
      fs.existsSync.mockReturnValue(false);

      loadEnvironmentVariables();

      expect(consoleWarnSpy).toHaveBeenCalledWith(
        '⚠️  .env file not found at: /path/to/.env'
      );
      expect(consoleWarnSpy).toHaveBeenCalledWith(
        expect.stringContaining('Create .env file with required variables')
      );
      expect(dotenvMock.config).not.toHaveBeenCalled();
    });

    it('should handle dotenv config errors', () => {
      servicesModule.detectCloudEnvironment.mockReturnValue(false);
      servicesModule.findEnvFilePath.mockReturnValue('/path/to/.env');
      fs.existsSync.mockReturnValue(true);
      dotenvMock.config.mockReturnValue({
        error: new Error('Failed to parse .env file')
      });

      loadEnvironmentVariables();

      expect(consoleErrorSpy).toHaveBeenCalledWith(
        '❌ Error loading .env file: Failed to parse .env file'
      );
      expect(servicesModule.validateCriticalEnvironmentVariables).not.toHaveBeenCalled();
    });
  });

  describe('Integration: Environment Loading Priority', () => {
    let dotenvMock;

    beforeEach(() => {
      dotenvMock = require('dotenv');
      dotenvMock.config.mockImplementation(() => {
        // Simulate dotenv setting environment variables
        process.env.FROM_DOTENV = 'dotenv-value';
        process.env.OVERRIDE_TEST = 'dotenv-override';
        return {};
      });
    });

    it('should prioritize shell environment variables over .env file', () => {
      // Set up shell environment variables first
      process.env.OVERRIDE_TEST = 'shell-value';
      process.env.SHELL_ONLY = 'shell-only-value';

      // Mock local development
      servicesModule.detectCloudEnvironment = jest.fn().mockReturnValue(false);
      servicesModule.findEnvFilePath = jest.fn().mockReturnValue('/path/to/.env');
      servicesModule.validateCriticalEnvironmentVariables = jest.fn();
      fs.existsSync.mockReturnValue(true);

      // Reset the global flag and re-run the loading logic
      delete global.__QUEB_ENV_LOADED__;
      require('../../src/core/services');

      // Shell environment variables should be preserved
      expect(process.env.SHELL_ONLY).toBe('shell-only-value');

      // .env should not override shell variables (override: false)
      expect(process.env.OVERRIDE_TEST).toBe('shell-value');

      // .env should add new variables
      expect(process.env.FROM_DOTENV).toBe('dotenv-value');
    });
  });
});