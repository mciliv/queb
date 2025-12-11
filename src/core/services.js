const fs = require('fs');
const path = require('path');
const dotenv = require('dotenv');

// ===== Environment Loading =====
// CRITICAL: Always load .env for local development
// Production deployments should use proper env vars, but local dev MUST have .env
if (!global.__QUEB_ENV_LOADED__) {
  // Detect if running in cloud/production environment
  const isCloud = !!(
    process.env.GAE_APPLICATION ||           // Google App Engine
    process.env.GOOGLE_CLOUD_PROJECT ||      // GCP
    process.env.K_SERVICE ||                 // Cloud Run
    process.env.FUNCTION_NAME ||             // Cloud Functions
    process.env.AWS_EXECUTION_ENV ||         // AWS Lambda
    process.env.VERCEL ||                    // Vercel
    process.env.NETLIFY ||                   // Netlify
    process.env.NODE_ENV === 'production'    // Explicit prod
  );

  // Load .env for local development
  if (!isCloud) {
    const envPath = path.resolve(process.cwd(), '.env');

    if (fs.existsSync(envPath)) {
      dotenv.config({ path: envPath, override: false });

      // Warn if required env vars missing after load
      if (!process.env.OPENAI_API_KEY) {
        console.warn('⚠️  OPENAI_API_KEY not found in .env');
      }
      if (!process.env.OPENAI_MODEL) {
        console.warn('⚠️  OPENAI_MODEL not found in .env');
      }
    } else {
      console.warn('⚠️  .env file not found at:', envPath);
      console.warn('   Create .env file with:');
      console.warn('   OPENAI_API_KEY=sk-...');
      console.warn('   OPENAI_MODEL=gpt-X');
    }
  }

  global.__QUEB_ENV_LOADED__ = true;
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

  // Do NOT change the model
  openai: {
    apiKey: getString('OPENAI_API_KEY'),
    model: getString('OPENAI_MODEL', 'gpt-5'),
    timeout: getNumber('OPENAI_TIMEOUT', 30000)
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
    if (!config.openai.apiKey) {
      return null;
    }
    return createOpenAIClient();
  }, {
    tags: ['ai', 'external']
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
      aiClient: await c.get('openaiClient'),
      molecularProcessor: await c.get('molecularProcessor'),
      nameResolver: await c.get('nameResolver'),
      promptEngine: await c.get('promptEngine'),
      errorHandler: await c.get('errorHandler'),
      logger: await c.get('logger'),
      config: await c.get('config')
    });
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
  createTestContainer
};