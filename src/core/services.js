const fs = require('fs');
const path = require('path');
const dotenv = require('dotenv');
const { log, warn, error } = require('./logger');

// ===== Environment Loading =====
// CRITICAL: Load environment variables from SHELL ENVIRONMENT first, then .env file for local development
// Shell environment variables ALWAYS take priority over .env file variables
// Production deployments should use proper env vars, but local dev can use .env as fallback
if (!global.__QUEB_ENV_LOADED__) {
  loadEnvironmentVariables();
  global.__QUEB_ENV_LOADED__ = true;
}

/**
 * Load environment variables from shell environment and .env file for local development.
 * This function is extracted for testability and clarity.
 *
 * Environment Variable Loading Strategy:
 * 1. Always prioritize shell environment variables (from system/shell)
 * 2. For local development: load .env file as fallback (does NOT override shell vars)
 * 3. For production/cloud: rely on shell environment variables only
 * 4. Warn about missing critical environment variables
 */
function loadEnvironmentVariables() {
  // Detect if running in cloud/production environment
  const isCloudEnvironment = detectCloudEnvironment();

  if (isCloudEnvironment) {
    log('🏭 Running in production/cloud environment - using SHELL ENVIRONMENT variables only');
    return;
  }

  // Load .env for local development (shell environment variables take priority)
  log('🏠 Running in local development - SHELL ENVIRONMENT takes priority, .env file provides fallbacks');
  const envFilePath = findEnvFilePath();

  if (!fs.existsSync(envFilePath)) {
    warn('⚠️  .env file not found', { envFilePath });
    warn('   Create .env file with required variables:');
    warn('   OPENAI_API_KEY=sk-...');
    warn('   OPENAI_MODEL=gpt-4');
    warn('   See .env.example for full list');
    return;
  }

  log('📄 Loading environment variables', { envFilePath });
  const result = dotenv.config({ path: envFilePath, override: false });

  if (result.error) {
    error('❌ Error loading .env file', result.error);
    return;
  }

  // Load hotel admin key if file exists
  loadHotelAdminKey();

  // Validate critical environment variables after loading
  validateCriticalEnvironmentVariables();
}

/**
 * Detect if the application is running in a cloud/production environment.
 * @returns {boolean} true if running in cloud/production, false for local development
 */
function detectCloudEnvironment() {
  return !!(
    process.env.GAE_APPLICATION ||           // Google App Engine
    process.env.GOOGLE_CLOUD_PROJECT ||      // GCP
    process.env.K_SERVICE ||                 // Cloud Run
    process.env.FUNCTION_NAME ||             // Cloud Functions
    process.env.AWS_EXECUTION_ENV ||         // AWS Lambda
    process.env.VERCEL ||                    // Vercel
    process.env.NETLIFY ||                   // Netlify
    process.env.NODE_ENV === 'production'    // Explicit production
  );
}

/**
 * Find the path to the .env file by locating the project root.
 * @returns {string} Absolute path to .env file
 */
function findEnvFilePath() {
  const projectRoot = findProjectRoot();
  return path.resolve(projectRoot, '.env');
}

/**
 * Load hotel admin key from file into environment variable.
 */
function loadHotelAdminKey() {
  const projectRoot = findProjectRoot();
  const keyFilePath = path.resolve(projectRoot, '.hotel_admin_key');

  if (fs.existsSync(keyFilePath)) {
    try {
      const keyContent = fs.readFileSync(keyFilePath, 'utf8').trim();
      // Skip comment lines and extract the key
      const lines = keyContent.split('\n').filter(line => line.trim() && !line.trim().startsWith('//'));
      const key = lines[0]?.trim();

      if (key) {
        process.env.HOTEL_ADMIN_KEY = key;
        log('🔑 Loaded hotel admin key from file');
      } else {
        warn('⚠️  Hotel admin key file exists but contains no valid key');
      }
    } catch (error) {
      error('❌ Error reading hotel admin key file', error);
    }
  }
}

/**
 * Find the project root directory (directory containing package.json).
 * @returns {string} Absolute path to project root
 */
function findProjectRoot() {
  let currentDir = __dirname;
  const root = path.parse(__dirname).root;

  while (currentDir !== root) {
    if (fs.existsSync(path.join(currentDir, 'package.json'))) {
      return currentDir;
    }
    currentDir = path.dirname(currentDir);
  }

  // Fallback: assume we're in src/core and go up two levels
  return path.resolve(__dirname, '../..');
}

/**
 * Validate that critical environment variables are present after loading.
 */
function validateCriticalEnvironmentVariables() {
  const criticalVars = [
    { key: 'OPENAI_API_KEY', description: 'OpenAI API key for AI services' },
    { key: 'OPENAI_MODEL', description: 'OpenAI model name (e.g., gpt-4)' }
  ];

  let missingVars = [];

  for (const { key, description } of criticalVars) {
    if (!process.env[key]) {
      missingVars.push({ key, description });
    }
  }

  if (missingVars.length > 0) {
    warn('⚠️  Missing critical environment variables in .env', { missingVars });
    for (const { key, description } of missingVars) {
      warn(`   - ${key}: ${description}`);
    }
    warn('   Add these to your .env file');
  } else {
    log('✅ All critical environment variables loaded successfully');
  }
}

// ===== Helpers =====
function getString(key, defaultValue = null) {
  return process.env[key] || defaultValue;
}

function getNumber(key, defaultValue = null) {
  const value = process.env[key];
  if (value === undefined || value === null) return defaultValue;
  const parsed = Number(value);
  return isNaN(parsed) ? defaultValue : parsed;
}

function getBoolean(key, defaultValue = false) {
  const value = process.env[key];
  if (value === undefined || value === null) return defaultValue;
  return value.toLowerCase() === 'true' || value === '1';
}

// ===== Configuration =====
const configData = {
  nodeEnv: getString('NODE_ENV', 'development'),
  port: getNumber('PORT', 8080),

  // AI Service Configuration - supports multiple providers
  ai: {
    provider: getString('AI_PROVIDER', 'openai'), // 'openai' or 'xai'
    openai: {
      apiKey: getString('OPENAI_API_KEY'),
      model: getString('OPENAI_MODEL', process.env.NODE_ENV === 'test' ? 'mock-model' : 'gpt-3.5-turbo'), // Use mock for tests
      timeout: getNumber('OPENAI_TIMEOUT', process.env.NODE_ENV === 'test' ? 1000 : 30000) // Fast timeout for tests
    },
    xai: {
      apiKey: getString('XAI_API_KEY'),
      model: getString('XAI_MODEL', 'grok-beta'),
      timeout: getNumber('XAI_TIMEOUT', process.env.NODE_ENV === 'test' ? 1000 : 30000) // Fast timeout for tests
    }
  },

  database: {
    enabled: getBoolean('DB_ENABLED', false),
    host: getString('DB_HOST', 'localhost'),
    port: getNumber('DB_PORT', 5432),
    name: getString('DB_NAME', 'mol_users'),
    user: getString('DB_USER', 'mol_user'),
    password: getString('DB_PASSWORD'),
    maxConnections: getNumber('DB_MAX_CONNECTIONS', 20),
    idleTimeout: getNumber('DB_IDLE_TIMEOUT', 30000),
    connectionTimeout: getNumber('DB_CONNECTION_TIMEOUT', 2000)
  },

  chem: {
    primary: getString('CHEM_PRIMARY', 'pubchem'),
    enableAlternates: getBoolean('CHEM_ENABLE_ALTERNATES', false)
  },

  cloud: {
    isCloudFunction: process.env.FUNCTION_NAME || process.env.FUNCTION_TARGET || process.env.K_SERVICE,
    isAppEngine: process.env.GAE_APPLICATION || process.env.GOOGLE_CLOUD_PROJECT,
    project: getString('GOOGLE_CLOUD_PROJECT'),
    region: getString('REGION', 'us-central1')
  },

  ssl: {
    certPath: getString('SSL_CERT_PATH'),
    keyPath: getString('SSL_KEY_PATH'),
    httpsPort: getNumber('HTTPS_PORT', 3001)
  },

  development: {
    isTest: process.env.NODE_ENV === 'test' || !!process.env.JEST_WORKER_ID,
    isIntegrationTest: getBoolean('INTEGRATION_TEST', false),
    enableDebug: getBoolean('DEBUG', false),
    enableScreenshots: getBoolean('ENABLE_ERROR_SCREENSHOTS', false)
  }
};

// Config object with .get() method for compatibility
const config = {
  ...configData,

  get(path = null) {
    if (!path) return configData;
    return path.split('.').reduce((obj, key) => obj?.[key], configData);
  },

  isProduction() {
    return configData.nodeEnv === 'production';
  },

  isDevelopment() {
    return configData.nodeEnv === 'development';
  },

  isTest() {
    return configData.development.isTest;
  },

  isServerless() {
    return configData.cloud.isCloudFunction || configData.cloud.isAppEngine;
  },

  getDatabaseConfig() {
    return configData.database;
  },

  validate() {
    const errors = [];

    if (this.isProduction()) {
      if (!configData.openai.apiKey) {
        errors.push('OPENAI_API_KEY is required in production');
      }

      if (configData.database.enabled && !configData.database.password) {
        errors.push('DB_PASSWORD is required when database is enabled in production');
      }
    }

    if (configData.port < 1 || configData.port > 65535) {
      errors.push(`Invalid PORT: ${configData.port}. Must be between 1 and 65535`);
    }

    if (errors.length > 0) {
      throw new Error('Configuration validation failed:\n' + errors.map(e => `  - ${e}`).join('\n'));
    }

    return true;
  },

  getDebugInfo() {
    return {
      nodeEnv: configData.nodeEnv,
      port: configData.port,
      hasOpenAIKey: !!configData.openai.apiKey,
      databaseEnabled: configData.database.enabled,
      isServerless: this.isServerless()
    };
  }
};

// ===== OpenAI Client Factory =====
function createOpenAIClient() {
  const OpenAI = require('openai').default || require('openai');
  if (!OpenAI) return null;

  return new OpenAI({
    apiKey: config.openai.apiKey,
    timeout: config.openai.timeout,
    maxRetries: 2
  });
}

// ===== Dependency Injection Container =====
const ServiceContainer = require('./ServiceContainer');
const PromptEngine = require('./PromptEngine');
const ErrorHandler = require('./ErrorHandler');

function createContainer(overrides = {}) {
  const MolecularPredictionService = require('./MolecularPredictionService');
  const container = new ServiceContainer();

  // ===== Core Services =====
  container.register('config', () => config, {
    tags: ['core']
  });

  container.register('logger', async () => {
    const FileLogger = require('../server/services/file-logger');
    return new FileLogger({
      logDir: config.get('logging.directory') || 'logs'
    });
  }, {
    tags: ['core', 'logging']
  });

  container.register('errorHandler', async (c) => {
    const handler = new ErrorHandler();
    handler.initialize(await c.get('logger'));
    return handler;
  }, {
    tags: ['core']
  });

  container.register('promptEngine', () => PromptEngine, {
    tags: ['core', 'ai']
  });

  container.register('openaiClient', async () => {
    if (!config.ai.openai.apiKey) {
      return null;
    }
    return createOpenAIClient();
  }, {
    tags: ['ai', 'external']
  });

  container.register('aiService', async () => {
    const AIService = require('../server/services/agents/AIService');

    return new AIService(); // Always loads from environment variables
  }, {
    tags: ['ai', 'service']
  });


  container.register('database', async () => {
    const Database = require('../server/services/database');

    if (!config.database.enabled) {
      return null;
    }

    return new Database(config.getDatabaseConfig());
  }, {
    tags: ['data', 'storage']
  });

  container.register('userService', async (c) => {
    const UserService = require('../server/services/user-service');
    const db = await c.get('database');

    if (!db) {
      const SimpleUserService = require('../server/services/simple-user-service');
      return new SimpleUserService();
    }

    return new UserService(db);
  }, {
    tags: ['data', 'business']
  });

  // ===== Chemistry Services =====
  container.register('nameResolver', () => {
    const resolver = require('../server/services/name-resolver');
    return resolver;
  }, {
    tags: ['chemistry', 'external']
  });

  container.register('molecularProcessor', async () => {
    const MolecularProcessor = require('../server/services/molecular-processor');
    return new MolecularProcessor(
      config.get('paths.sdfDirectory')
    );
  }, {
    tags: ['chemistry', 'processing']
  });

  // ===== Business Services =====
  container.register('structuralizer', async (c) => {
    const Structuralizer = require('../server/services/Structuralizer');

    return new Structuralizer({
      aiService: await c.get('aiService'),
      molecularProcessor: await c.get('molecularProcessor'),
      nameResolver: await c.get('nameResolver'),
      promptEngine: await c.get('promptEngine'),
      errorHandler: await c.get('errorHandler'),
      logger: await c.get('logger'),
      config: await c.get('config')
    });
  });

  // ===== Logging Services =====
  // PORT: Abstract interface for logging (domain depends on this)
  container.register('logErrorPort', async (c) => {
    const LogErrorPort = require('../server/services/ports/log-error-port');
    // PORT CONTRACT: Domain only knows the interface, not implementation
    return LogErrorPort; // Return the class, not instance (ports are abstract)
  }, {
    tags: ['ports', 'logging']
  });

  // ADAPTER: Concrete implementation of LogErrorPort using FileLogger
  container.register('fileLoggerAdapter', async (c) => {
    const FileLoggerAdapter = require('../server/services/adapters/file-logger-adapter');
    const fileLogger = await c.get('logger'); // Inject the FileLogger instance

    // ADAPTER PATTERN: Wraps infrastructure (FileLogger) to match domain contract (LogErrorPort)
    return new FileLoggerAdapter(fileLogger);
  }, {
    tags: ['adapters', 'logging', 'infrastructure']
  });

  // USE-CASE: Domain logic orchestrator that depends on ports
  container.register('logErrorUseCase', async (c) => {
    const LogErrorUseCase = require('../server/services/log-error-use-case');

    // DEPENDENCY INJECTION: Use-case gets port implementation (could be swapped)
    const logErrorPort = await c.get('fileLoggerAdapter'); // Could be any LogErrorPort implementation

    // USE-CASE PATTERN: Domain logic with injected infrastructure dependencies
    return new LogErrorUseCase({ logErrorPort });
  }, {
    tags: ['use-cases', 'logging', 'domain']
  });

  container.register('molecularPredictionService', async (c) => {
    const service = new MolecularPredictionService({
      aiClient: await c.get('openaiClient'),
      molecularProcessor: await c.get('molecularProcessor'),
      nameResolver: await c.get('nameResolver'),
      logger: await c.get('logger')
    });

    await service.initialize();
    return service;
  }, {
    tags: ['business', 'prediction']
  });

  // ===== Performance Services =====
  container.register('cache', async () => {
    const PerformanceCache = require('../server/services/performance-cache');

    return new PerformanceCache({
      maxSize: config.get('cache.maxSize') || 100,
      ttl: config.get('cache.ttl') || 300000
    });
  }, {
    tags: ['performance']
  });

  // ===== Apply overrides =====
  Object.entries(overrides).forEach(([name, service]) => {
    container.register(name, service);
  });

  return container;
}

function createTestContainer(mocks = {}) {
  const defaults = {
    logger: {
      info: jest.fn(),
      error: jest.fn(),
      warn: jest.fn(),
      debug: jest.fn()
    },
    config: {
      get: jest.fn((key) => {
        const testConfig = {
          'openai.apiKey': 'test-key',
          'database.enabled': false,
          'paths.sdfDirectory': '/tmp/test-sdf'
        };
        return testConfig[key];
      }),
      getDatabaseConfig: jest.fn(),
      isProduction: jest.fn(() => false),
      isDevelopment: jest.fn(() => true)
    },
    openaiClient: {
      chat: {
        completions: {
          create: jest.fn()
        }
      }
    }
  };

  return createContainer({ ...defaults, ...mocks });
}

module.exports = {
  config,
  createOpenAIClient,
  createContainer,
  createTestContainer,
  // Environment loading functions (exported for testing)
  loadEnvironmentVariables,
  detectCloudEnvironment,
  findEnvFilePath,
  findProjectRoot,
  validateCriticalEnvironmentVariables
};