#!/usr/bin/env node

/**
 * PostgreSQL Database Restore - Node.js version
 * Replaces database/scripts/restore.sh
 */

const { execSync } = require('child_process');
const fs = require('fs');
const path = require('path');
const readline = require('readline');

// Colored logging utilities (replacing bash utils.sh)
const log = {
  info: (msg) => console.log(`\x1b[32m✅ ${msg}\x1b[0m`),
  warn: (msg) => console.log(`\x1b[33m⚠️  ${msg}\x1b[0m`),
  error: (msg) => console.log(`\x1b[31m❌ ${msg}\x1b[0m`)
};

function getDbConfig() {
  return {
    host: process.env.DB_HOST || 'localhost',
    port: parseInt(process.env.DB_PORT, 10) || 5432,
    database: process.env.DB_NAME || 'mol_users',
    user: process.env.DB_USER || 'mol_user',
    password: process.env.DB_PASSWORD
  };
}

function showUsage() {
  log.error('Usage: node restore.js <backup_file.sql>');
  log.info('Available backups:');
  
  try {
    const backupDir = './backups';
    if (fs.existsSync(backupDir)) {
      const files = fs.readdirSync(backupDir)
        .filter(file => file.startsWith('backup_') && file.endsWith('.sql'))
        .map(file => {
          const filePath = path.join(backupDir, file);
          const stats = fs.statSync(filePath);
          return {
            name: file,
            path: filePath,
            modified: stats.mtime.toISOString().replace('T', ' ').split('.')[0]
          };
        })
        .sort((a, b) => b.modified.localeCompare(a.modified));

      if (files.length > 0) {
        files.forEach(file => {
          console.log(`  ${file.path} (${file.modified})`);
        });
      } else {
        log.info('No backups found in ./backups/');
      }
    } else {
      log.info('No backups directory found');
    }
  } catch (error) {
    log.warn(`Could not list backups: ${error.message}`);
  }
}

async function askConfirmation(question) {
  // In test environment, auto-confirm to avoid hanging
  if (process.env.NODE_ENV === 'test') {
    return true;
  }

  const rl = readline.createInterface({
    input: process.stdin,
    output: process.stdout
  });

  return new Promise((resolve) => {
    rl.question(question, (answer) => {
      rl.close();
      resolve(answer.toLowerCase() === 'y' || answer.toLowerCase() === 'yes');
    });
  });
}

async function restoreDatabase(backupFile) {
  // Check if backup file exists
  if (!fs.existsSync(backupFile)) {
    log.error(`Backup file not found: ${backupFile}`);
    process.exit(1);
  }

  log.warn('⚠️  This will overwrite the current database!');
  const confirmed = await askConfirmation(`Are you sure you want to restore from ${backupFile}? (y/N): `);
  
  if (!confirmed) {
    log.info('Restore cancelled');
    process.exit(0);
  }

  log.info(`Starting PostgreSQL restore from: ${backupFile}`);

  const dbConfig = getDbConfig();
  const env = { 
    ...process.env, 
    PGPASSWORD: dbConfig.password 
  };

  try {
    // Terminate active connections to the database
    log.info('Terminating active connections...');
    try {
      const terminateCommand = `psql -h ${dbConfig.host} -p ${dbConfig.port} -U ${dbConfig.user} -d postgres -c "SELECT pg_terminate_backend(pid) FROM pg_stat_activity WHERE datname = '${dbConfig.database}' AND pid <> pg_backend_pid();"`;
      execSync(terminateCommand, { env, stdio: 'pipe' });
    } catch (error) {
      log.warn('Could not terminate connections (may not be necessary)');
    }

    // Drop and recreate database
    log.info('Recreating database...');
    const recreateCommand = `psql -h ${dbConfig.host} -p ${dbConfig.port} -U ${dbConfig.user} -d postgres -c "DROP DATABASE IF EXISTS ${dbConfig.database}; CREATE DATABASE ${dbConfig.database} OWNER ${dbConfig.user};"`;
    
    try {
      execSync(recreateCommand, { env, stdio: 'pipe' });
    } catch (error) {
      log.error('Failed to recreate database');
      throw error;
    }

    // Restore from backup
    log.info('Restoring from backup...');
    const restoreCommand = `psql -h ${dbConfig.host} -p ${dbConfig.port} -U ${dbConfig.user} -d ${dbConfig.database} -f "${backupFile}"`;
    
    execSync(restoreCommand, { 
      env,
      stdio: 'inherit'
    });

    log.info('✅ Database restore completed successfully');

  } catch (error) {
    log.error(`Restore failed: ${error.message}`);
    process.exit(1);
  }
}

async function main() {
  const args = process.argv.slice(2);
  
  if (args.length === 0) {
    showUsage();
    process.exit(1);
  }

  const backupFile = args[0];
  await restoreDatabase(backupFile);
}

if (require.main === module) {
  main().catch(error => {
    log.error(`Restore process failed: ${error.message}`);
    process.exit(1);
  });
}

module.exports = { restoreDatabase, log };
