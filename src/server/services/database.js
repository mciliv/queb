// Database configuration and initialization
const { services } = require('../../core/services');

let pool = null;
let dbConnected = false;

const initializeDatabase = () => {
  const dbConfig = config.getDatabaseConfig();
  if (!dbConfig.enabled) {
    console.log('💡 Database disabled - running without user storage');
    return null;
  }

  try {
    const { Pool } = require('pg');

    const poolConfig = {
      host: dbConfig.host,
      port: dbConfig.port,
      database: dbConfig.name,
      user: dbConfig.user,
      password: dbConfig.password,
      max: dbConfig.maxConnections,
      idleTimeoutMillis: dbConfig.idleTimeout,
      connectionTimeoutMillis: dbConfig.connectionTimeout,
    };

    pool = new Pool(poolConfig);

    // Database connection error handling
    pool.on('error', (err, client) => {
      console.error('🔴 Unexpected error on idle client', err);
      console.log('💡 Database connection will be retried automatically');
    });

    console.log('✅ PostgreSQL module loaded successfully');
    return pool;
  } catch (error) {
    console.warn('⚠️ PostgreSQL module not available - running without database');
    console.log('💡 Install with: npm install pg pg-pool');
    return null;
  }
};

const testDatabaseConnection = async () => {
  if (!pool) {
    console.warn('⚠️ Database not available - running without persistent user storage');
    return false;
  }

  try {
    const client = await pool.connect();
    const result = await client.query('SELECT NOW()');
    client.release();
    console.log('✅ Database connected successfully');
    dbConnected = true;
    return true;
  } catch (err) {
    console.error('🔴 Database connection failed:', err.message);
    dbConnected = false;
    return false;
  }
};

const setupDatabase = async (UserService) => {
  const dbPool = initializeDatabase();
  const userService = (dbPool && UserService) ? new UserService(dbPool) : null;

  if (!userService) {
    console.log('⚠️ User service not available - running without persistent user storage');
    return { pool: null, userService: null };
  }

  try {
    const connected = await testDatabaseConnection();
    if (connected) {
      await userService.initializeTables();
    } else {
      console.log('⚠️ Database not connected - running without persistent user storage');
    }
  } catch (error) {
    console.error('🔴 Database initialization failed:', error.message);
    console.log('💡 Server will continue but user data will not persist');
  }

  return { pool: dbPool, userService };
};

module.exports = {
  initializeDatabase,
  testDatabaseConnection,
  setupDatabase,
  getPool: () => pool,
  isConnected: () => dbConnected,
};