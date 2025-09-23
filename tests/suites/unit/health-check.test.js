/**
 * Unit tests for database health-check.js script
 */

const { checkDatabaseHealth, log } = require('../../../database/scripts/health-check.js');

// Mock pg Client
const mockClient = {
  connect: jest.fn(),
  query: jest.fn(),
  end: jest.fn()
};

jest.mock('pg', () => ({
  Client: jest.fn(() => mockClient)
}));

describe('health-check.js', () => {
  let consoleSpy;
  let processExitSpy;

  beforeEach(() => {
    consoleSpy = jest.spyOn(console, 'log').mockImplementation();
    processExitSpy = jest.spyOn(process, 'exit').mockImplementation();
    
    // Reset all mocks
    mockClient.connect.mockClear();
    mockClient.query.mockClear();
    mockClient.end.mockClear();
  });

  afterEach(() => {
    consoleSpy.mockRestore();
    processExitSpy.mockRestore();
  });

  describe('logging utilities', () => {
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

  describe('checkDatabaseHealth function', () => {
    test('should successfully complete health check with all queries', async () => {
      // Mock successful connection and queries
      mockClient.connect.mockResolvedValueOnce();
      mockClient.query
        .mockResolvedValueOnce() // Initial connection test
        .mockResolvedValueOnce({ rows: [{ size: '10 MB' }] }) // Database size
        .mockResolvedValueOnce({ rows: [{ count: '5' }] }) // Active connections
        .mockResolvedValueOnce({ rows: [{ count: '3' }] }) // Table count
        .mockResolvedValueOnce({ rows: [{ count: '2' }] }); // Recent activity
      mockClient.end.mockResolvedValueOnce();

      await checkDatabaseHealth();

      expect(mockClient.connect).toHaveBeenCalledTimes(1);
      expect(mockClient.query).toHaveBeenCalledTimes(5);
      expect(mockClient.end).toHaveBeenCalledTimes(1);
      expect(processExitSpy).not.toHaveBeenCalled();

      // Verify log messages
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ 🔍 Checking PostgreSQL database health...\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ ✅ Database connection: OK\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ 📊 Database size: 10 MB\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ 🔗 Active connections: 5\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ 📋 Tables: 3\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ 📈 Recent activity (24h): 2 users\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ 🎉 Health check completed!\x1b[0m');
    });

    test('should handle connection failure', async () => {
      const connectionError = new Error('Connection failed');
      mockClient.connect.mockRejectedValueOnce(connectionError);
      mockClient.end.mockResolvedValueOnce();

      await checkDatabaseHealth();

      expect(mockClient.connect).toHaveBeenCalledTimes(1);
      expect(mockClient.end).toHaveBeenCalledTimes(1);
      expect(processExitSpy).toHaveBeenCalledWith(1);
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[31m❌ ❌ Database connection: FAILED - Connection failed\x1b[0m');
    });

    test('should handle individual query failures gracefully', async () => {
      // Mock successful connection but failed queries
      mockClient.connect.mockResolvedValueOnce();
      mockClient.query
        .mockResolvedValueOnce() // Initial connection test
        .mockRejectedValueOnce(new Error('Size query failed')) // Database size fails
        .mockResolvedValueOnce({ rows: [{ count: '5' }] }) // Active connections succeeds
        .mockRejectedValueOnce(new Error('Table query failed')) // Table count fails
        .mockRejectedValueOnce(new Error('Activity query failed')); // Recent activity fails
      mockClient.end.mockResolvedValueOnce();

      await checkDatabaseHealth();

      expect(mockClient.connect).toHaveBeenCalledTimes(1);
      expect(mockClient.query).toHaveBeenCalledTimes(5);
      expect(mockClient.end).toHaveBeenCalledTimes(1);
      expect(processExitSpy).not.toHaveBeenCalled();

      // Should show warnings for failed queries
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[33m⚠️  ⚠️  Could not determine database size\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ 🔗 Active connections: 5\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[33m⚠️  ⚠️  Could not count tables\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ 📈 Recent activity: Unable to check (users table may not exist)\x1b[0m');
    });

    test('should handle empty query results', async () => {
      // Mock successful connection but empty results
      mockClient.connect.mockResolvedValueOnce();
      mockClient.query
        .mockResolvedValueOnce() // Initial connection test
        .mockResolvedValueOnce({ rows: [] }) // Database size - empty
        .mockResolvedValueOnce({ rows: [] }) // Active connections - empty
        .mockResolvedValueOnce({ rows: [] }) // Table count - empty
        .mockResolvedValueOnce({ rows: [] }); // Recent activity - empty
      mockClient.end.mockResolvedValueOnce();

      await checkDatabaseHealth();

      expect(processExitSpy).not.toHaveBeenCalled();
      // Empty results don't generate warnings in the actual implementation
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ ✅ Database connection: OK\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ 🎉 Health check completed!\x1b[0m');
    });
  });

  describe('database configuration', () => {
    const originalEnv = process.env;

    beforeEach(() => {
      process.env = { ...originalEnv };
    });

    afterEach(() => {
      process.env = originalEnv;
    });

    test('should use environment variables when available', () => {
      process.env.DB_HOST = 'test-host';
      process.env.DB_PORT = '5433';
      process.env.DB_NAME = 'test_db';
      process.env.DB_USER = 'test_user';
      process.env.DB_PASSWORD = 'test_pass';

      // Re-require the module to pick up new env vars
      delete require.cache[require.resolve('../../../database/scripts/health-check.js')];
      const { checkDatabaseHealth } = require('../../../database/scripts/health-check.js');

      // The config should be used in the Client constructor
      expect(require('pg').Client).toHaveBeenCalledWith({
        host: 'test-host',
        port: 5433,
        database: 'test_db',
        user: 'test_user',
        password: 'test_pass'
      });
    });

    test('should use default values when environment variables are not set', () => {
      delete process.env.DB_HOST;
      delete process.env.DB_PORT;
      delete process.env.DB_NAME;
      delete process.env.DB_USER;
      delete process.env.DB_PASSWORD;

      // Re-require the module to pick up new env vars
      delete require.cache[require.resolve('../../../database/scripts/health-check.js')];
      const { checkDatabaseHealth } = require('../../../database/scripts/health-check.js');

      expect(require('pg').Client).toHaveBeenCalledWith({
        host: 'localhost',
        port: 5432,
        database: 'mol_users',
        user: 'mol_user',
        password: undefined
      });
    });
  });
});
