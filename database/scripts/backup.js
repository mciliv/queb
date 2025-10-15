#!/usr/bin/env node

/**
 * PostgreSQL Database Backup - Node.js version  
 * Replaces database/scripts/backup.sh
 */

const { execSync } = require('child_process');
const fs = require('fs');
const path = require('path');

// Colored logging utilities (replacing bash utils.sh)
const log = {
  info: (msg) => console.log(`\x1b[32m✅ ${msg}\x1b[0m`),
  warn: (msg) => console.log(`\x1b[33m⚠️  ${msg}\x1b[0m`),
  error: (msg) => console.log(`\x1b[31m❌ ${msg}\x1b[0m`)
};

<<<<<<< Updated upstream
function getDbConfig() {
  return {
    host: process.env.DB_HOST || 'localhost',
    port: process.env.DB_PORT || 5432,
=======
// Database configuration - read dynamically to support test environment
function getDbConfig() {
  return {
    host: process.env.DB_HOST || 'localhost',
    port: parseInt(process.env.DB_PORT) || 5432,
>>>>>>> Stashed changes
    database: process.env.DB_NAME || 'mol_users',
    user: process.env.DB_USER || 'mol_user',
    password: process.env.DB_PASSWORD
  };
}

<<<<<<< Updated upstream
function getBackupDir() {
  return process.env.BACKUP_DIR || './backups';
}

function getBackupFilePath(now = new Date()) {
  const backupDir = getBackupDir();
  const timestamp = now.toISOString()
    .replace(/[:.]/g, '-')
    .replace('T', '_')
    .split('.')[0];
  return path.join(backupDir, `backup_${timestamp}.sql`);
=======
// Backup configuration - read dynamically
function getBackupConfig() {
  const backupDir = process.env.BACKUP_DIR || './backups';
  const timestamp = new Date().toISOString()
    .replace(/[:.]/g, '-')
    .replace('T', '_')
    .split('.')[0]; // Format: YYYY-MM-DD_HH-MM-SS

  const backupFile = path.join(backupDir, `backup_${timestamp}.sql`);
  
  return { backupDir, backupFile, timestamp };
>>>>>>> Stashed changes
}

async function createBackup() {
  const dbConfig = getDbConfig();
  const { backupDir, backupFile } = getBackupConfig();
  log.info('Starting PostgreSQL backup...');

  try {
    const dbConfig = getDbConfig();
    const backupDir = getBackupDir();
    const backupFile = getBackupFilePath();
    // Create backup directory if it doesn't exist
    if (!fs.existsSync(backupDir)) {
      fs.mkdirSync(backupDir, { recursive: true });
    }

    // Perform backup using pg_dump
    log.info(`Creating backup: ${backupFile}`);
    
    const pgDumpCommand = [
      'pg_dump',
      `-h ${dbConfig.host}`,
      `-p ${dbConfig.port}`,
      `-U ${dbConfig.user}`,
      `-d ${dbConfig.database}`,
      `-f "${backupFile}"`
    ].join(' ');

    // Set password environment variable for pg_dump
    const env = { 
      ...process.env, 
      PGPASSWORD: dbConfig.password 
    };

    execSync(pgDumpCommand, { 
      env,
      stdio: 'inherit'
    });

    log.info(`✅ Backup completed: ${backupFile}`);

    // Clean up old backups (keep last 7 days)
    await cleanupOldBackups();

  } catch (error) {
    log.error(`Backup failed: ${error.message}`);
    process.exit(1);
  }
}

async function cleanupOldBackups() {
  const { backupDir } = getBackupConfig();
  try {
    const backupDir = getBackupDir();
    const files = fs.readdirSync(backupDir);
    const backupFiles = files.filter(file => file.startsWith('backup_') && file.endsWith('.sql'));
    
    const sevenDaysAgo = new Date();
    sevenDaysAgo.setDate(sevenDaysAgo.getDate() - 7);

    let cleanedCount = 0;
    
    for (const file of backupFiles) {
      const filePath = path.join(backupDir, file);
      const stats = fs.statSync(filePath);
      
      if (stats.mtime < sevenDaysAgo) {
        fs.unlinkSync(filePath);
        cleanedCount++;
      }
    }

    if (cleanedCount > 0) {
      log.info(`🧹 Cleaned up ${cleanedCount} old backup(s)`);
    }
  } catch (error) {
    log.warn(`Could not clean up old backups: ${error.message}`);
  }
}

if (require.main === module) {
  createBackup().catch(error => {
    log.error(`Backup process failed: ${error.message}`);
    process.exit(1);
  });
}

module.exports = { createBackup, cleanupOldBackups, log };
