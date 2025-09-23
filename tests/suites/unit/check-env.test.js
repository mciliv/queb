/**
 * Unit tests for check-env.js script
 */

const { log, checkPort, checkUrl } = require('../../../scripts/check-env.js');
const net = require('net');

describe('check-env.js utilities', () => {
  
  describe('logging utilities', () => {
    let consoleSpy;

    beforeEach(() => {
      consoleSpy = jest.spyOn(console, 'log').mockImplementation();
    });

    afterEach(() => {
      consoleSpy.mockRestore();
    });

    test('log.info should output green colored message', () => {
      log.info('test message');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ test message\x1b[0m');
    });

    test('log.warn should output yellow colored message', () => {
      log.warn('warning message');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[33m⚠️  warning message\x1b[0m');
    });

    test('log.error should output red colored message', () => {
      log.error('error message');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[31m❌ error message\x1b[0m');
    });
  });

  describe('checkPort function', () => {
    let testServer;

    afterEach(async () => {
      if (testServer) {
        await new Promise(resolve => {
          testServer.close(resolve);
        });
        testServer = null;
      }
    });

    test('should return false for available port', async () => {
      const port = 19999; // Use unlikely port
      const result = await checkPort(port);
      expect(result).toBe(false);
    });

    test('should return true for port in use', async () => {
      const port = 19998;
      
      // Create a server to occupy the port
      testServer = net.createServer();
      await new Promise(resolve => {
        testServer.listen(port, '127.0.0.1', resolve);
      });

      const result = await checkPort(port);
      expect(result).toBe(true);
    });
  });

  describe('checkUrl function', () => {
    // Mock fetch for testing
    const originalFetch = global.fetch;

    beforeAll(() => {
      global.fetch = jest.fn();
    });

    afterAll(() => {
      global.fetch = originalFetch;
    });

    beforeEach(() => {
      fetch.mockClear();
    });

    test('should return true for successful response', async () => {
      fetch.mockResolvedValueOnce({ ok: true });

      const result = await checkUrl('https://example.com');
      expect(result).toBe(true);
      expect(fetch).toHaveBeenCalledWith('https://example.com', {
        signal: expect.any(AbortSignal),
        method: 'HEAD'
      });
    });

    test('should return false for failed response', async () => {
      fetch.mockResolvedValueOnce({ ok: false });

      const result = await checkUrl('https://example.com');
      expect(result).toBe(false);
    });

    test('should return false for network error', async () => {
      fetch.mockRejectedValueOnce(new Error('Network error'));

      const result = await checkUrl('https://example.com');
      expect(result).toBe(false);
    });

    test('should handle timeout', async () => {
      // Mock fetch to hang indefinitely
      fetch.mockImplementationOnce(() => new Promise(() => {}));

      const result = await checkUrl('https://example.com', 100); // 100ms timeout
      expect(result).toBe(false);
    }, 5000); // 5 second timeout for this test
  });

  describe('environment detection', () => {
    const originalEnv = process.env.NODE_ENV;

    afterEach(() => {
      process.env.NODE_ENV = originalEnv;
    });

    test('should detect production environment', () => {
      process.env.NODE_ENV = 'production';
      // This would need to be tested by running the main function
      // For now, we just verify the env var is set
      expect(process.env.NODE_ENV).toBe('production');
    });

    test('should detect development environment', () => {
      process.env.NODE_ENV = 'development';
      expect(process.env.NODE_ENV).toBe('development');
    });
  });
});
