/**
 * Unit tests for database backup.js script
 */

const { createBackup, cleanupOldBackups, log } = require('../../../database/scripts/backup.js');
const { execSync } = require('child_process');
const fs = require('fs');
const path = require('path');

// Mock dependencies
jest.mock('child_process');
jest.mock('fs');
jest.mock('path');

describe('backup.js', () => {
  let consoleSpy;
  let processExitSpy;
  const originalEnv = process.env;

  beforeEach(() => {
    consoleSpy = jest.spyOn(console, 'log').mockImplementation();
    processExitSpy = jest.spyOn(process, 'exit').mockImplementation();
    
    // Reset all mocks
    execSync.mockClear();
    fs.existsSync.mockClear();
    fs.mkdirSync.mockClear();
    fs.readdirSync.mockClear();
    fs.statSync.mockClear();
    fs.unlinkSync.mockClear();
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

  describe('createBackup function', () => {
    beforeEach(() => {
      // Mock path.join to return predictable paths
      path.join.mockImplementation((...args) => args.join('/'));
    });

    test('should successfully create backup', async () => {
      // Mock successful backup creation
      fs.existsSync.mockReturnValueOnce(false); // Backup dir doesn't exist
      fs.mkdirSync.mockReturnValueOnce(); // Create directory succeeds
      execSync.mockReturnValueOnce(); // pg_dump succeeds
      fs.readdirSync.mockReturnValueOnce([]); // No old backups to clean

      await createBackup();

      expect(fs.mkdirSync).toHaveBeenCalledWith('./backups', { recursive: true });
      expect(execSync).toHaveBeenCalledWith(
        expect.stringContaining('pg_dump'),
        expect.objectContaining({
          env: expect.objectContaining({
            PGPASSWORD: undefined // Default password
          }),
          stdio: 'inherit'
        })
      );
      expect(processExitSpy).not.toHaveBeenCalled();
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ Starting PostgreSQL backup...\x1b[0m');
    });

    test('should use custom backup directory from environment', async () => {
      process.env.BACKUP_DIR = '/custom/backup/path';
      
      fs.existsSync.mockReturnValueOnce(true); // Directory exists
      execSync.mockReturnValueOnce(); // pg_dump succeeds
      fs.readdirSync.mockReturnValueOnce([]); // No old backups

      await createBackup();

      expect(fs.mkdirSync).not.toHaveBeenCalled(); // Directory already exists
      expect(path.join).toHaveBeenCalledWith('/custom/backup/path', expect.any(String));
    });

    test('should use database config from environment variables', async () => {
      process.env.DB_HOST = 'test-host';
      process.env.DB_PORT = '5433';
      process.env.DB_NAME = 'test_db';
      process.env.DB_USER = 'test_user';
      process.env.DB_PASSWORD = 'secret123';

      fs.existsSync.mockReturnValueOnce(true);
      execSync.mockReturnValueOnce();
      fs.readdirSync.mockReturnValueOnce([]);

      await createBackup();

      expect(execSync).toHaveBeenCalledWith(
        expect.stringContaining('-h test-host -p 5433 -U test_user -d test_db'),
        expect.objectContaining({
          env: expect.objectContaining({
            PGPASSWORD: 'secret123'
          }),
          stdio: 'inherit'
        })
      );
    });

    test('should handle backup failure', async () => {
      fs.existsSync.mockReturnValueOnce(true);
      execSync.mockImplementationOnce(() => {
        throw new Error('pg_dump failed');
      });

      await createBackup();

      expect(processExitSpy).toHaveBeenCalledWith(1);
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[31m❌ Backup failed: pg_dump failed\x1b[0m');
    });

    test('should create backup directory if it does not exist', async () => {
      fs.existsSync.mockReturnValueOnce(false); // Directory doesn't exist
      fs.mkdirSync.mockReturnValueOnce(); // Create succeeds
      execSync.mockReturnValueOnce(); // Backup succeeds
      fs.readdirSync.mockReturnValueOnce([]);

      await createBackup();

      expect(fs.mkdirSync).toHaveBeenCalledWith('./backups', { recursive: true });
    });

    test('should include timestamp in backup filename', async () => {
      fs.existsSync.mockReturnValueOnce(true);
      execSync.mockReturnValueOnce();
      fs.readdirSync.mockReturnValueOnce([]);

      // Mock Date to return predictable timestamp
      const mockDate = new Date('2023-01-15T10:30:45.123Z');
      jest.spyOn(global, 'Date').mockImplementation(() => mockDate);

      await createBackup();

      expect(execSync).toHaveBeenCalledWith(
        expect.stringContaining('backup_2023-01-15_10-30-45.sql'),
        expect.any(Object)
      );

      global.Date.mockRestore();
    });
  });

  describe('cleanupOldBackups function', () => {
    beforeEach(() => {
      path.join.mockImplementation((...args) => args.join('/'));
    });

    test('should remove old backup files', async () => {
      const oldDate = new Date();
      oldDate.setDate(oldDate.getDate() - 10); // 10 days old

      const recentDate = new Date();
      recentDate.setDate(recentDate.getDate() - 3); // 3 days old

      fs.readdirSync.mockReturnValueOnce([
        'backup_old.sql',
        'backup_recent.sql',
        'other_file.txt', // Should be ignored
        'backup_another_old.sql'
      ]);

      fs.statSync
        .mockReturnValueOnce({ mtime: oldDate }) // old backup
        .mockReturnValueOnce({ mtime: recentDate }) // recent backup
        .mockReturnValueOnce({ mtime: oldDate }); // another old backup

      fs.unlinkSync.mockReturnValue();

      await cleanupOldBackups();

      expect(fs.unlinkSync).toHaveBeenCalledTimes(2);
      expect(fs.unlinkSync).toHaveBeenCalledWith('./backups/backup_old.sql');
      expect(fs.unlinkSync).toHaveBeenCalledWith('./backups/backup_another_old.sql');
      expect(consoleSpy).toHaveBeenCalledWith('\x1b[32m✅ 🧹 Cleaned up 2 old backup(s)\x1b[0m');
    });

    test('should handle cleanup errors gracefully', async () => {
      fs.readdirSync.mockImplementationOnce(() => {
        throw new Error('Cannot read directory');
      });

      await cleanupOldBackups();

      expect(consoleSpy).toHaveBeenCalledWith('\x1b[33m⚠️  Could not clean up old backups: Cannot read directory\x1b[0m');
      expect(processExitSpy).not.toHaveBeenCalled();
    });

    test('should not remove recent backups', async () => {
      const recentDate = new Date();
      recentDate.setDate(recentDate.getDate() - 3); // 3 days old (within 7 days)

      fs.readdirSync.mockReturnValueOnce(['backup_recent.sql']);
      fs.statSync.mockReturnValueOnce({ mtime: recentDate });

      await cleanupOldBackups();

      expect(fs.unlinkSync).not.toHaveBeenCalled();
      expect(consoleSpy).not.toHaveBeenCalledWith(expect.stringContaining('Cleaned up'));
    });

    test('should only process backup files with correct naming pattern', async () => {
      const oldDate = new Date();
      oldDate.setDate(oldDate.getDate() - 10);

      fs.readdirSync.mockReturnValueOnce([
        'backup_old.sql', // Should be processed
        'other_backup.sql', // Should be ignored (wrong prefix)
        'backup_recent.txt', // Should be ignored (wrong extension)
        'random_file.sql' // Should be ignored (wrong prefix)
      ]);

      fs.statSync.mockReturnValueOnce({ mtime: oldDate });

      await cleanupOldBackups();

      expect(fs.statSync).toHaveBeenCalledTimes(1); // Only called for backup_old.sql
      expect(fs.unlinkSync).toHaveBeenCalledTimes(1);
    });
  });
});
