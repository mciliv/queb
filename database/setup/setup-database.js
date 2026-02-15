#!/usr/bin/env node

/**
 * PostgreSQL Database Setup - Node.js version
 * Replaces database/setup/setup-database.sh
 */

const { execSync, exec } = require('child_process');
const { Client } = require('pg');
const fs = require('fs');
const path = require('path');
const { promisify } = require('util');

const execAsync = promisify(exec);

// Colored logging utilities (replacing bash utils.sh)
const log = {
  info: (msg) => console.log(`\x1b[32m✅ ${msg}\x1b[0m`),
  warn: (msg) => console.log(`\x1b[33m⚠️  ${msg}\x1b[0m`),
  error: (msg) => console.log(`\x1b[31m❌ ${msg}\x1b[0m`)
};

// Default database configuration
const defaultConfig = {
  host: 'localhost',
  port: 5432,
  database: 'mol_users',
  user: 'mol_user',
  password: 'mol_password'
};

// Check if PostgreSQL is installed and running
async function checkPostgreSQL() {
  log.info('Checking PostgreSQL installation...');

  // Check if PostgreSQL client is installed
  try {
    execSync('psql --version', { stdio: 'pipe' });
    log.info('PostgreSQL client found');
  } catch (error) {
    log.error('PostgreSQL client not found');
    log.info('Please install PostgreSQL:');
    console.log('  macOS: brew install postgresql');
    console.log('  Ubuntu: sudo apt-get install postgresql postgresql-contrib');
    console.log('  Windows: Download from https://www.postgresql.org/download/');
    process.exit(1);
  }

  // Check if PostgreSQL server is running
  try {
    execSync(`pg_isready -h ${defaultConfig.host} -p ${defaultConfig.port}`, { stdio: 'pipe' });
    log.info('PostgreSQL server is running');
  } catch (error) {
    log.error('PostgreSQL server is not running');
    log.info('Start PostgreSQL server:');
    console.log('  macOS: brew services start postgresql');
    console.log('  Ubuntu: sudo systemctl start postgresql');
    console.log('  Windows: Start from Services or pgAdmin');
    process.exit(1);
  }
}

// Get current user for PostgreSQL connection
function getCurrentUser() {
  try {
    return execSync('whoami', { encoding: 'utf8' }).trim();
  } catch (error) {
    return 'postgres';
  }
}

// Create database and user
async function setupDatabase() {
  log.info('Setting up database and user...');

  const currentUser = getCurrentUser();
  let postgresUser = 'postgres';

  // Try to determine the correct superuser
  try {
    // Test if we can connect as postgres user
    const testClient = new Client({
      host: defaultConfig.host,
      port: defaultConfig.port,
      user: 'postgres',
      database: 'postgres'
    });
    await testClient.connect();
    await testClient.end();
    postgresUser = 'postgres';
  } catch (error) {
    // Try with current user (common on macOS with Homebrew)
    postgresUser = currentUser;
    log.info(`Attempting to connect as user: ${postgresUser}`);
  }

  const adminClient = new Client({
    host: defaultConfig.host,
    port: defaultConfig.port,
    user: postgresUser,
    database: 'postgres'
  });

  try {
    await adminClient.connect();

    // Create user if it doesn't exist
    log.info(`Creating database user: ${defaultConfig.user}`);
    await adminClient.query(`
      DO $$
      BEGIN
        IF NOT EXISTS (SELECT 1 FROM pg_roles WHERE rolname = $1) THEN
          CREATE USER ${defaultConfig.user} WITH PASSWORD $2;
          ALTER USER ${defaultConfig.user} CREATEDB;
        END IF;
      END
      $$;
    `, [defaultConfig.user, defaultConfig.password]);

    // Create database if it doesn't exist
    log.info(`Creating database: ${defaultConfig.database}`);
    const dbExistsResult = await adminClient.query(
      'SELECT 1 FROM pg_database WHERE datname = $1',
      [defaultConfig.database]
    );

    if (dbExistsResult.rows.length === 0) {
      await adminClient.query(`CREATE DATABASE ${defaultConfig.database} OWNER ${defaultConfig.user}`);
    } else {
      log.info('Database already exists');
    }

    // Grant privileges
    log.info('Granting privileges...');
    await adminClient.query(`GRANT ALL PRIVILEGES ON DATABASE ${defaultConfig.database} TO ${defaultConfig.user}`);

  } catch (error) {
    log.warn(`Could not create user/database with ${postgresUser}, trying alternative method...`);
    
    // Try creating database directly with current user
    try {
      execSync(`createdb ${defaultConfig.database}`, { stdio: 'pipe' });
      log.info('Database created with current user');
    } catch (createError) {
      log.error('Failed to create database. Please run manually:');
      console.log(`  sudo -u postgres createuser -s ${defaultConfig.user}`);
      console.log(`  sudo -u postgres createdb -O ${defaultConfig.user} ${defaultConfig.database}`);
      process.exit(1);
    }
  } finally {
    await adminClient.end();
  }
}

// Test database connection
async function testConnection() {
  log.info('Testing database connection...');

  const testClient = new Client({
    host: defaultConfig.host,
    port: defaultConfig.port,
    user: defaultConfig.user,
    password: defaultConfig.password,
    database: defaultConfig.database
  });

  try {
    await testClient.connect();
    await testClient.query('SELECT NOW()');
    log.info('Database connection successful');
  } catch (error) {
    log.error('Database connection failed');
    log.info('Please check your PostgreSQL installation and try again');
    process.exit(1);
  } finally {
    await testClient.end();
  }
}

function createEnvFile() {
  log.info('Creating .env file for local development...');

  // Find queb project root (directory containing package.json)
  const findProjectRoot = () => {
    let currentDir = __dirname;
    while (currentDir !== path.dirname(currentDir)) {
      if (fs.existsSync(path.join(currentDir, 'package.json'))) {
        return currentDir;
      }
      currentDir = path.dirname(currentDir);
    }
    return path.resolve(__dirname, '../..'); // Fallback: go up from database/setup to queb
  };
  const projectRoot = findProjectRoot();
  const envFile = path.join(projectRoot, '.env');

  if (!fs.existsSync(envFile)) {
    const envContent = `# PostgreSQL Configuration for Local Development
DB_HOST=${defaultConfig.host}
DB_PORT=${defaultConfig.port}
DB_NAME=${defaultConfig.database}
DB_USER=${defaultConfig.user}
DB_PASSWORD=${defaultConfig.password}
`;

    fs.writeFileSync(envFile, envContent);
    log.info(`.env file created at ${envFile}`);
    log.warn('Please update the API keys in .env file');
  } else {
    log.info('.env file already exists, skipping creation');
  }
}

// Main setup function
async function main() {
  console.log('🧬 Molecular Analysis App - Database Setup');
  console.log('==========================================');

  try {
    await checkPostgreSQL();
    await setupDatabase();
    await testConnection();
    createEnvFile();

    console.log('');
    log.info('Database setup completed successfully!');
    console.log('');
    log.info('Connection details:');
    console.log(`  Host: ${defaultConfig.host}`);
    console.log(`  Port: ${defaultConfig.port}`);
    console.log(`  Database: ${defaultConfig.database}`);
    console.log(`  User: ${defaultConfig.user}`);
    console.log('');
    log.info('To connect manually:');
    console.log(`  PGPASSWORD=${defaultConfig.password} psql -h ${defaultConfig.host} -p ${defaultConfig.port} -U ${defaultConfig.user} -d ${defaultConfig.database}`);
    console.log('');
    log.info('Your app will automatically create the required tables on startup.');

  } catch (error) {
    log.error(`Setup failed: ${error.message}`);
    process.exit(1);
  }
}

if (require.main === module) {
  main();
}

module.exports = { 
  checkPostgreSQL, 
  setupDatabase, 
  testConnection, 
  createEnvFile, 
  log,
  defaultConfig 
};
