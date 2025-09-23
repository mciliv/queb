/**
 * Unit tests for database restore.js script
 */

const { restoreDatabase, log } = require('../../../database/scripts/restore.js');
const { execSync } = require('child_process');
const fs = require('fs');
const readline = require('readline');

// Mock dependencies
jest.mock('child_process');
jest.mock('fs');
jest.mock('readline');

describe('restore.js', () => {
  let consoleSpy;
  let processExitSpy;
  const originalEnv = process.env;

  beforeEach(() => {
    consoleSpy = jest.spyOn(console, 'log').mockImplementation();
    processExitSpy = jest.spyOn(process, 'exit').mockImplementation();
    
    // Reset all mocks
    execSync.mockClear();
    fs.existsSync.mockClear();
    fs.readdirSync.mockClear();
    fs.statSync.mockClear();
    readline.createInterface.mockClear();

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

  describe('restoreDatabase function', () => {
    const mockRl = {
      question: jest.fn(),
      close: jest.fn()
    };

    beforeEach(() => {
      readline.createInterface.mockReturnValue(mockRl);
      mockRl.question.mockClear();
      mockRl.close.mockClear();
    });

    test('should exit if backup file does not exist', async () => {
      fs.existsSync.mockReturnValueOnce(false);

      await restoreDatabase('/path/to/nonexistent.sql');

      expect(processExitSpy).toHaveBeenCalledWith(1);
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[31m❌ Backup file not found: /path/to/nonexistent.sql\x1b[0m');
    });

    test('should exit if user cancels confirmation', async () => {
      fs.existsSync.mockReturnValueOnce(true);
      
      // Mock user saying 'no'
      mockRl.question.mockImplementationOnce((question, callback) => {
        callback('n');
      });

      await restoreDatabase('/path/to/backup.sql');

      expect(processExitSpy).toHaveBeenCalledWith(0);
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ Restore cancelled\x1b[0m');
    });

    test('should successfully restore database when user confirms', async () => {
      fs.existsSync.mockReturnValueOnce(true);
      
      // Mock user saying 'yes'
      mockRl.question.mockImplementationOnce((question, callback) => {
        callback('y');
      });

      // Mock successful execSync calls
      execSync
        .mockReturnValueOnce() // Terminate connections
        .mockReturnValueOnce() // Drop/recreate database
        .mockReturnValueOnce(); // Restore from backup

      await restoreDatabase('/path/to/backup.sql');

      expect(execSync).toHaveBeenCalledTimes(3);
      expect(processExitSpy).not.toHaveBeenCalled();
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ ✅ Database restore completed successfully\x1b[0m');
    });

    test('should handle different confirmation responses', async () => {
      const testCases = [
        { input: 'Y', shouldProceed: true },
        { input: 'yes', shouldProceed: true },
        { input: 'YES', shouldProceed: true },
        { input: 'n', shouldProceed: false },
        { input: 'no', shouldProceed: false },
        { input: 'maybe', shouldProceed: false },
        { input: '', shouldProceed: false }
      ];

      for (const testCase of testCases) {
        // Reset mocks
        processExitSpy.mockClear();
        consoleSpy.mockClear();
        execSync.mockClear();

        fs.existsSync.mockReturnValueOnce(true);
        mockRl.question.mockImplementationOnce((question, callback) => {
          callback(testCase.input);
        });

        if (testCase.shouldProceed) {
          execSync
            .mockReturnValueOnce()
            .mockReturnValueOnce()
            .mockReturnValueOnce();
        }

        await restoreDatabase('/path/to/backup.sql');

        if (testCase.shouldProceed) {
          expect(execSync).toHaveBeenCalledTimes(3);
          expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ ✅ Database restore completed successfully\x1b[0m');
        } else {
          expect(processExitSpy).toHaveBeenCalledWith(0);
          expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ Restore cancelled\x1b[0m');
        }
      }
    });

    test('should use database config from environment variables', async () => {
      process.env.DB_HOST = 'test-host';
      process.env.DB_PORT = '5433';
      process.env.DB_NAME = 'test_db';
      process.env.DB_USER = 'test_user';
      process.env.DB_PASSWORD = 'secret123';

      fs.existsSync.mockReturnValueOnce(true);
      mockRl.question.mockImplementationOnce((question, callback) => {
        callback('y');
      });

      execSync
        .mockReturnValueOnce()
        .mockReturnValueOnce()
        .mockReturnValueOnce();

      await restoreDatabase('/path/to/backup.sql');

      // Check that the correct database config was used in execSync calls
      expect(execSync).toHaveBeenCalledWith(
        expect.stringContaining('-h test-host -p 5433 -U test_user -d postgres'),
        expect.objectContaining({
          env: expect.objectContaining({
            PGPASSWORD: 'secret123'
          })
        })
      );

      expect(execSync).toHaveBeenCalledWith(
        expect.stringContaining('CREATE DATABASE test_db OWNER test_user'),
        expect.any(Object)
      );

      expect(execSync).toHaveBeenCalledWith(
        expect.stringContaining('-h test-host -p 5433 -U test_user -d test_db'),
        expect.any(Object)
      );
    });

    test('should handle connection termination failure gracefully', async () => {
      fs.existsSync.mockReturnValueOnce(true);
      mockRl.question.mockImplementationOnce((question, callback) => {
        callback('y');
      });

      execSync
        .mockImplementationOnce(() => { throw new Error('Connection termination failed'); })
        .mockReturnValueOnce() // Drop/recreate database succeeds
        .mockReturnValueOnce(); // Restore succeeds

      await restoreDatabase('/path/to/backup.sql');

      expect(consoleSpy).toHaveBeenCalledWith('\x1b[33m⚠️  Could not terminate connections (may not be necessary)\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ ✅ Database restore completed successfully\x1b[0m');
      expect(processExitSpy).not.toHaveBeenCalled();
    });

    test('should handle database recreation failure', async () => {
      fs.existsSync.mockReturnValueOnce(true);
      mockRl.question.mockImplementationOnce((question, callback) => {
        callback('y');
      });

      execSync
        .mockReturnValueOnce() // Connection termination succeeds
        .mockImplementationOnce(() => { throw new Error('Database recreation failed'); });

      await restoreDatabase('/path/to/backup.sql');

      expect(consoleSpy).toHaveBeenCalledWith('\x1b[31m❌ Failed to recreate database\x1b[0m');
      expect(processExitSpy).toHaveBeenCalledWith(1);
    });

    test('should handle restore failure', async () => {
      fs.existsSync.mockReturnValueOnce(true);
      mockRl.question.mockImplementationOnce((question, callback) => {
        callback('y');
      });

      execSync
        .mockReturnValueOnce() // Connection termination succeeds
        .mockReturnValueOnce() // Database recreation succeeds
        .mockImplementationOnce(() => { throw new Error('Restore failed'); });

      await restoreDatabase('/path/to/backup.sql');

      expect(consoleSpy).toHaveBeenCalledWith('\x1b[31m❌ Restore failed: Restore failed\x1b[0m');
      expect(processExitSpy).toHaveBeenCalledWith(1);
    });

    test('should show warning about overwriting database', async () => {
      fs.existsSync.mockReturnValueOnce(true);
      mockRl.question.mockImplementationOnce((question, callback) => {
        expect(question).toContain('This will overwrite the current database!');
        expect(question).toContain('/path/to/backup.sql');
        callback('y');
      });

      execSync
        .mockReturnValueOnce()
        .mockReturnValueOnce()
        .mockReturnValueOnce();

      await restoreDatabase('/path/to/backup.sql');

      expect(consoleSpy).toHaveBeenCalledWith('\x1b[33m⚠️  ⚠️  This will overwrite the current database!\x1b[0m');
    });

    test('should log all restoration steps', async () => {
      fs.existsSync.mockReturnValueOnce(true);
      mockRl.question.mockImplementationOnce((question, callback) => {
        callback('y');
      });

      execSync
        .mockReturnValueOnce()
        .mockReturnValueOnce()
        .mockReturnValueOnce();

      await restoreDatabase('/path/to/backup.sql');

      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ Starting PostgreSQL restore from: /path/to/backup.sql\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ Terminating active connections...\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ Recreating database...\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ Restoring from backup...\x1b[0m');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ ✅ Database restore completed successfully\x1b[0m');
    });
  });
});
