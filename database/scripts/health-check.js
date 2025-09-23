#!/usr/bin/env node

/**
 * PostgreSQL Database Health Check - Node.js version
 * Replaces database/scripts/health-check.sh
 */

const { Client } = require('pg');

// Colored logging utilities (replacing bash utils.sh)
const log = {
  info: (msg) => console.log(`\x1b[32m✅ ${msg}\x1b[0m`),
  warn: (msg) => console.log(`\x1b[33m⚠️  ${msg}\x1b[0m`),
  error: (msg) => console.log(`\x1b[31m❌ ${msg}\x1b[0m`)
};

// Database configuration
const dbConfig = {
  host: process.env.DB_HOST || 'localhost',
  port: process.env.DB_PORT || 5432,
  database: process.env.DB_NAME || 'mol_users',
  user: process.env.DB_USER || 'mol_user',
  password: process.env.DB_PASSWORD
};

async function checkDatabaseHealth() {
  log.info('🔍 Checking PostgreSQL database health...');

  const client = new Client(dbConfig);

  try {
    // Check connection
    await client.connect();
    await client.query('SELECT 1');
    log.info('✅ Database connection: OK');

    // Check database size
    try {
      const sizeResult = await client.query(`
        SELECT pg_size_pretty(pg_database_size($1)) as size
      `, [dbConfig.database]);
      
      if (sizeResult.rows.length > 0) {
        log.info(`📊 Database size: ${sizeResult.rows[0].size}`);
      }
    } catch (error) {
      log.warn('⚠️  Could not determine database size');
    }

    // Check active connections
    try {
      const connectionsResult = await client.query(`
        SELECT COUNT(*) as count FROM pg_stat_activity WHERE datname = $1
      `, [dbConfig.database]);
      
      if (connectionsResult.rows.length > 0) {
        log.info(`🔗 Active connections: ${connectionsResult.rows[0].count}`);
      }
    } catch (error) {
      log.warn('⚠️  Could not check connections');
    }

    // Check table count
    try {
      const tablesResult = await client.query(`
        SELECT COUNT(*) as count FROM information_schema.tables 
        WHERE table_schema = 'public'
      `);
      
      if (tablesResult.rows.length > 0) {
        log.info(`📋 Tables: ${tablesResult.rows[0].count}`);
      }
    } catch (error) {
      log.warn('⚠️  Could not count tables');
    }

    // Check recent activity (last 24 hours)
    try {
      const activityResult = await client.query(`
        SELECT COUNT(*) as count FROM users 
        WHERE last_used > NOW() - INTERVAL '24 hours'
      `);
      
      if (activityResult.rows.length > 0) {
        log.info(`📈 Recent activity (24h): ${activityResult.rows[0].count} users`);
      }
    } catch (error) {
      log.info('📈 Recent activity: Unable to check (users table may not exist)');
    }

    log.info('🎉 Health check completed!');

  } catch (error) {
    log.error(`❌ Database connection: FAILED - ${error.message}`);
    process.exit(1);
  } finally {
    await client.end();
  }
}

if (require.main === module) {
  checkDatabaseHealth().catch(error => {
    log.error(`Health check failed: ${error.message}`);
    process.exit(1);
  });
}

module.exports = { checkDatabaseHealth, log };
