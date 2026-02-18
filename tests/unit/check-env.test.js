const net = require('net');

// Define a persistent mock server object
const mockServer = {
  once: jest.fn(),
  listen: jest.fn(),
  close: jest.fn(),
};

// Mock the 'net' module to always return the same mockServer instance
jest.mock('net', () => ({
  createServer: jest.fn(() => mockServer),
}));

/**
 * checkPort function from the original check-env.js
 */
async function checkPort(port) {
  return new Promise((resolve) => {
    const server = net.createServer();
    
    server.once('error', (err) => {
      if (err.code === 'EADDRINUSE') {
        resolve(true); // Port is in use
      } else {
        resolve(false);
      }
    });
    
    server.once('listening', () => {
      server.close();
      resolve(false); // Port is available
    });
    
    server.listen(port, '127.0.0.1');
  });
}

/**
 * checkUrl function from the original check-env.js
 */
async function checkUrl(url, timeout = 5000) {
  const controller = new AbortController();
  let timeoutId = null;
  try {
    timeoutId = setTimeout(() => controller.abort(), timeout);

    const response = await fetch(url, {
      signal: controller.signal,
      method: 'HEAD'
    });

    return response.ok;
  } catch (error) {
    return false;
  } finally {
    if (timeoutId) clearTimeout(timeoutId);
  }
}

describe('Environment Utility Checks', () => {
  beforeEach(() => {
    // Reset mocks before each test
    jest.clearAllMocks();

    // Reset the implementation of the mock server's methods
    mockServer.once.mockImplementation((event, callback) => {
      mockServer.once[event] = callback;
    });
    mockServer.listen.mockImplementation((port, host) => {
      if (mockServer.listen.shouldError) {
        if (mockServer.once.error) {
          mockServer.once.error({ code: 'EADDRINUSE' });
        }
      } else {
        if (mockServer.once.listening) {
          mockServer.once.listening();
        }
      }
    });
  });

  describe('checkPort', () => {
    test('should return true if port is in use', async () => {
      mockServer.listen.shouldError = true;
      const result = await checkPort(3000);
      expect(result).toBe(true);
      expect(net.createServer).toHaveBeenCalledTimes(1);
      expect(mockServer.listen).toHaveBeenCalledWith(3000, '127.0.0.1');
      expect(mockServer.close).not.toHaveBeenCalled();
    });

    test('should return false if port is available', async () => {
      mockServer.listen.shouldError = false;
      const result = await checkPort(3000);
      expect(result).toBe(false);
      expect(net.createServer).toHaveBeenCalledTimes(1);
      expect(mockServer.listen).toHaveBeenCalledWith(3000, '127.0.0.1');
      expect(mockServer.close).toHaveBeenCalledTimes(1);
    });
  });

  describe('checkUrl', () => {
    test('should return true if URL is accessible (response ok)', async () => {
      global.fetch.mockResolvedValueOnce({ ok: true });
      const result = await checkUrl('http://example.com');
      expect(result).toBe(true);
      expect(global.fetch).toHaveBeenCalledTimes(1);
      expect(global.fetch).toHaveBeenCalledWith('http://example.com', expect.any(Object));
    });

    test('should return false if URL is not accessible (response not ok)', async () => {
      global.fetch.mockResolvedValueOnce({ ok: false });
      const result = await checkUrl('http://example.com');
      expect(result).toBe(false);
      expect(global.fetch).toHaveBeenCalledTimes(1);
    });

    test('should return false if fetch throws an error', async () => {
      global.fetch.mockRejectedValueOnce(new Error('Network error'));
      const result = await checkUrl('http://example.com');
      expect(result).toBe(false);
      expect(global.fetch).toHaveBeenCalledTimes(1);
    });

    test('should return false if fetch times out', async () => {
      jest.useFakeTimers(); // Mock timers
      global.fetch.mockImplementationOnce(({ signal }) => {
        return new Promise((resolve, reject) => {
          if (signal) {
            signal.addEventListener('abort', () => {
              reject(new DOMException('Aborted', 'AbortError'));
            });
          }
        });
      });
      
      const checkUrlPromise = checkUrl('http://example.com', 10); // Start checkUrl
      jest.runAllTimers(); // Advance all timers, triggering controller.abort()
      const result = await checkUrlPromise; // Await the result

      expect(result).toBe(false);
      expect(global.fetch).toHaveBeenCalledTimes(1);
      jest.useRealTimers(); // Restore real timers
    });
  });
});
