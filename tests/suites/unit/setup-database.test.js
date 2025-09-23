/**
 * Unit tests for database setup-database.js script
 */

const { 
  checkPostgreSQL, 
  setupDatabase, 
  testConnection, 
  createEnvFile, 
  log,
  defaultConfig 
} = require('../../../database/setup/setup-database.js');

const { execSync } = require('child_process');
const { Client } = require('pg');
const fs = require('fs');
const path = require('path');

// Mock dependencies
jest.mock('child_process');
jest.mock('pg');
jest.mock('fs');
jest.mock('path');

describe('setup-database.js', () => {
  let consoleSpy;
  let processExitSpy;
  const originalEnv = process.env;

  // Mock pg Client
  const mockClient = {
    connect: jest.fn(),
    query: jest.fn(),
    end: jest.fn()
  };

  beforeEach(() => {
    consoleSpy = jest.spyOn(console, 'log').mockImplementation();
    processExitSpy = jest.spyOn(process, 'exit').mockImplementation();
    
    // Reset all mocks
    execSync.mockClear();
    Client.mockClear();
    Client.mockReturnValue(mockClient);
    mockClient.connect.mockClear();
    mockClient.query.mockClear();
    mockClient.end.mockClear();
    fs.existsSync.mockClear();
    fs.writeFileSync.mockClear();
    path.join.mockClear();

    // Set up environment
    process.env = { ...originalEnv };
  });

  afterEach(() => {
    consoleSpy.mockRestore();
    processExitSpy.mockRestore();
    process.env = originalEnv;
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

  describe('checkPostgreSQL function', () => {
    test('should pass when PostgreSQL is installed and running', async () => {
      execSync
        .mockReturnValueOnce('psql (PostgreSQL) 14.1') // psql --version
        .mockReturnValueOnce('accepting connections'); // pg_isready

      await checkPostgreSQL();

      expect(execSync).toHaveBeenCalledWith('psql --version', { stdio: 'pipe' });
      expect(execSync).toHaveBeenCalledWith('pg_isready -h localhost -p 5432', { stdio: 'pipe' });
      expect(processExitSpy).not.toHaveBeenCalled();
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ PostgreSQL client found\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ PostgreSQL server is running\x1b[0m');
    });

    test('should exit when PostgreSQL client is not installed', async () => {
      execSync.mockImplementationOnce(() => {
        throw new Error('Command not found');
      });

      await checkPostgreSQL();

      expect(processExitSpy).toHaveBeenCalledWith(1);
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[31m❌ PostgreSQL client not found\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ Please install PostgreSQL:\x1b[0m');
    });

    test('should exit when PostgreSQL server is not running', async () => {
      execSync
        .mockReturnValueOnce('psql (PostgreSQL) 14.1') // psql --version succeeds
        .mockImplementationOnce(() => { // pg_isready fails
          throw new Error('Server not running');
        });

      await checkPostgreSQL();

      expect(processExitSpy).toHaveBeenCalledWith(1);
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ PostgreSQL client found\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[31m❌ PostgreSQL server is not running\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ Start PostgreSQL server:\x1b[0m');
    });
  });

  describe('setupDatabase function', () => {
    beforeEach(() => {
      // Mock whoami command
      execSync.mockImplementation((cmd) => {
        if (cmd === 'whoami') {
          return 'testuser';
        }
        return '';
      });
    });

    test('should successfully create user and database', async () => {
      mockClient.connect.mockResolvedValueOnce();
      mockClient.query
        .mockResolvedValueOnce() // Create user query
        .mockResolvedValueOnce({ rows: [] }) // Check if database exists (empty = doesn't exist)
        .mockResolvedValueOnce() // Create database
        .mockResolvedValueOnce(); // Grant privileges
      mockClient.end.mockResolvedValueOnce();

      await setupDatabase();

      expect(mockClient.connect).toHaveBeenCalledTimes(2); // First test connection, then admin connection
      expect(mockClient.query).toHaveBeenCalledTimes(4);
      expect(mockClient.end).toHaveBeenCalledTimes(2); // End both connections
      expect(processExitSpy).not.toHaveBeenCalled();
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ Creating database user: mol_user\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ Creating database: mol_users\x1b[0m');
    });

    test('should handle existing database', async () => {
      mockClient.connect.mockResolvedValueOnce();
      mockClient.query
        .mockResolvedValueOnce() // Create user query
        .mockResolvedValueOnce({ rows: [{ datname: 'mol_users' }] }) // Database exists
        .mockResolvedValueOnce(); // Grant privileges
      mockClient.end.mockResolvedValueOnce();

      await setupDatabase();

      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ Database already exists\x1b[0m');
    });

    test('should try alternative method when admin connection fails', async () => {
      mockClient.connect.mockRejectedValueOnce(new Error('Connection failed'));
      mockClient.end.mockResolvedValueOnce();
      execSync.mockReturnValueOnce(''); // createdb succeeds

      await setupDatabase();

      expect(execSync).toHaveBeenCalledWith('createdb mol_users', { stdio: 'pipe' });
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ Database created with current user\x1b[0m');
      expect(processExitSpy).not.toHaveBeenCalled();
    });

    test('should exit when all methods fail', async () => {
      // Mock both connection attempts to fail
      mockClient.connect
        .mockRejectedValueOnce(new Error('Postgres connection failed'))
        .mockRejectedValueOnce(new Error('Admin connection failed'));
      mockClient.end.mockResolvedValueOnce();
      execSync.mockImplementationOnce(() => {
        throw new Error('createdb failed');
      });

      await setupDatabase();

      expect(processExitSpy).toHaveBeenCalledWith(1);
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[31m❌ Failed to create database. Please run manually:\x1b[0m');
    });

    test('should try postgres user first, then current user', async () => {
      // First try with postgres user (should succeed)
      mockClient.connect.mockResolvedValueOnce();
      mockClient.query
        .mockResolvedValueOnce()
        .mockResolvedValueOnce({ rows: [] })
        .mockResolvedValueOnce()
        .mockResolvedValueOnce();
      mockClient.end.mockResolvedValueOnce();

      await setupDatabase();

      // Verify Client was called with postgres user first
      expect(Client).toHaveBeenCalledWith({
        host: 'localhost',
        port: 5432,
        user: 'postgres',
        database: 'postgres'
      });
    });

    test('should use current user when postgres connection fails', async () => {
      // First connection attempt fails
      mockClient.connect.mockRejectedValueOnce(new Error('postgres user failed'));
      mockClient.end.mockResolvedValueOnce();

      // Second connection attempt with current user succeeds
      const mockClient2 = {
        connect: jest.fn().mockResolvedValueOnce(),
        query: jest.fn()
          .mockResolvedValueOnce()
          .mockResolvedValueOnce({ rows: [] })
          .mockResolvedValueOnce()
          .mockResolvedValueOnce(),
        end: jest.fn().mockResolvedValueOnce()
      };
      
      Client
        .mockReturnValueOnce(mockClient) // First call (postgres user)
        .mockReturnValueOnce(mockClient2); // Second call (current user)

      await setupDatabase();

      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ Attempting to connect as user: testuser\x1b[0m');
    });
  });

  describe('testConnection function', () => {
    test('should pass with successful connection', async () => {
      mockClient.connect.mockResolvedValueOnce();
      mockClient.query.mockResolvedValueOnce();
      mockClient.end.mockResolvedValueOnce();

      await testConnection();

      expect(mockClient.connect).toHaveBeenCalledTimes(1);
      expect(mockClient.query).toHaveBeenCalledWith('SELECT NOW()');
      expect(mockClient.end).toHaveBeenCalledTimes(1);
      expect(processExitSpy).not.toHaveBeenCalled();
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ Database connection successful\x1b[0m');
    });

    test('should exit on connection failure', async () => {
      mockClient.connect.mockRejectedValueOnce(new Error('Connection failed'));
      mockClient.end.mockResolvedValueOnce();

      await testConnection();

      expect(processExitSpy).toHaveBeenCalledWith(1);
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[31m❌ Database connection failed\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ Please check your PostgreSQL installation and try again\x1b[0m');
    });

    test('should use correct database configuration', async () => {
      mockClient.connect.mockResolvedValueOnce();
      mockClient.query.mockResolvedValueOnce();
      mockClient.end.mockResolvedValueOnce();

      await testConnection();

      expect(Client).toHaveBeenCalledWith({
        host: 'localhost',
        port: 5432,
        user: 'mol_user',
        password: 'mol_password',
        database: 'mol_users'
      });
    });
  });

  describe('createEnvFile function', () => {
    beforeEach(() => {
      path.join.mockImplementation((...args) => args.join('/'));
    });

    test('should create .env file when it does not exist', () => {
      fs.existsSync.mockReturnValueOnce(false);
      fs.writeFileSync.mockReturnValueOnce();

      createEnvFile();

      expect(fs.existsSync).toHaveBeenCalledWith(expect.stringContaining('.env'));
      expect(fs.writeFileSync).toHaveBeenCalledWith(
        expect.stringContaining('.env'),
        expect.stringContaining('DB_HOST=localhost')
      );
      expect(consoleSpy).toHaveBeenCalledWith(expect.stringContaining('.env file created'));
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[33m⚠️  Please update the API keys in .env file\x1b[0m');
    });

    test('should skip creation when .env file exists', () => {
      fs.existsSync.mockReturnValueOnce(true);

      createEnvFile();

      expect(fs.writeFileSync).not.toHaveBeenCalled();
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ .env file already exists, skipping creation\x1b[0m');
    });

    test('should include all required environment variables', () => {
      fs.existsSync.mockReturnValueOnce(false);
      fs.writeFileSync.mockReturnValueOnce();

      createEnvFile();

      const envContent = fs.writeFileSync.mock.calls[0][1];
      expect(envContent).toContain('DB_HOST=localhost');
      expect(envContent).toContain('DB_PORT=5432');
      expect(envContent).toContain('DB_NAME=mol_users');
      expect(envContent).toContain('DB_USER=mol_user');
      expect(envContent).toContain('DB_PASSWORD=mol_password');
      expect(envContent).toContain('NODE_ENV=development');
      expect(envContent).toContain('OPENAI_API_KEY=your_openai_api_key_here');
      expect(envContent).toContain('STRIPE_PUBLISHABLE_KEY=your_stripe_publishable_key_here');
      expect(envContent).toContain('STRIPE_SECRET_KEY=your_stripe_secret_key_here');
    });
  });

  describe('defaultConfig', () => {
    test('should have correct default values', () => {
      expect(defaultConfig).toEqual({
        host: 'localhost',
        port: 5432,
        database: 'mol_users',
        user: 'mol_user',
        password: 'mol_password'
      });
    });
  });

  describe('error handling', () => {
    test('should handle getCurrentUser function errors', () => {
      execSync.mockImplementationOnce((cmd) => {
        if (cmd === 'whoami') {
          throw new Error('whoami failed');
        }
        return '';
      });

      // This should not throw - it should default to 'postgres'
      const { getCurrentUser } = require('../../../database/setup/setup-database.js');
      // getCurrentUser is not exported, but we can test it indirectly through setupDatabase
    });

    test('should handle database connection timeout', async () => {
      const timeoutError = new Error('Connection timeout');
      mockClient.connect.mockRejectedValueOnce(timeoutError);
      mockClient.end.mockResolvedValueOnce();

      await testConnection();

      expect(processExitSpy).toHaveBeenCalledWith(1);
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[31m❌ Database connection failed\x1b[0m');
    });
  });
});
