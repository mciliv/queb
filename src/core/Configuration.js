/**
 * Unified Configuration System
 * 
 * This module applies the "Deep Modules" principle by providing a simple interface
 * that hides the complexity of environment variable management, validation, and defaults.
 * 
 * Philosophy: "The best modules are those whose interfaces are much simpler than 
 * their implementations" - John Ousterhout
 */

const path = require('path');

class Configuration {
  constructor() {
    this._config = null;
    this._initialized = false;
    this._loadEnvironment();
  }

  /**
   * Load environment variables with fallback hierarchy
   * This method hides the complexity of multiple .env file loading
   */
  _loadEnvironment() {
    try {
      // Load dotenv with multiple file support
      require('dotenv').config({ path: '.env.defaults', override: false });
      require('dotenv').config({ path: '.env', override: true });
      require('dotenv').config({ path: '.env.local', override: true });
      
      this._initialized = true;
    } catch (error) {
      console.warn('⚠️ Environment loading failed, using defaults:', error.message);
      this._initialized = true;
    }
  }

  /**
   * Get configuration with type safety and validation
   * Deep module: complex validation logic hidden behind simple interface
   */
  _buildConfig() {
    if (this._config) return this._config;

    this._config = {
      // Environment
      nodeEnv: this._getString('NODE_ENV', 'development'),
      port: this._getNumber('PORT', 8080),
      
      // AI Services
      openai: {
        apiKey: this._getString('OPENAI_API_KEY'),
        model: this._getString('OPENAI_MODEL', 'gpt-4o-mini'),
        timeout: this._getNumber('OPENAI_TIMEOUT', 30000)
      },

      // Database
      database: {
        enabled: this._getBoolean('DB_ENABLED', false),
        host: this._getString('DB_HOST', 'localhost'),
        port: this._getNumber('DB_PORT', 5432),
        name: this._getString('DB_NAME', 'mol_users'),
        user: this._getString('DB_USER', 'mol_user'),
        password: this._getString('DB_PASSWORD'),
        maxConnections: this._getNumber('DB_MAX_CONNECTIONS', 20),
        idleTimeout: this._getNumber('DB_IDLE_TIMEOUT', 30000),
        connectionTimeout: this._getNumber('DB_CONNECTION_TIMEOUT', 2000)
      },

      // Payments
      payments: {
        enabled: this._getBoolean('PAYMENTS_ENABLED', false),
        devMode: this._getBoolean('PAYMENTS_DEV_MODE', true),
        required: this._getBoolean('PAYMENTS_REQUIRED', false),
        stripePublishableKey: this._getString('STRIPE_PUBLISHABLE_KEY', 'pk_test_demo_key_for_development')
      },

      // Chemistry providers
      chem: {
        primary: this._getString('CHEM_PRIMARY', 'pubchem'), // 'pubchem' | 'chembl'
        enableAlternates: this._getBoolean('CHEM_ENABLE_ALTERNATES', false)
      },

      // Cloud Services
      cloud: {
        isCloudFunction: !!(process.env.FUNCTION_NAME || process.env.FUNCTION_TARGET || process.env.K_SERVICE),
        isAppEngine: !!(process.env.GAE_APPLICATION || process.env.GOOGLE_CLOUD_PROJECT),
        isNetlify: !!process.env.NETLIFY,
        project: this._getString('GOOGLE_CLOUD_PROJECT'),
        region: this._getString('REGION', 'us-central1')
      },

      // SSL
      ssl: {
        certPath: this._getString('SSL_CERT_PATH'),
        keyPath: this._getString('SSL_KEY_PATH'),
        httpsPort: this._getNumber('HTTPS_PORT', 3001)
      },

      // Development
      development: {
        isTest: process.env.NODE_ENV === 'test' || !!process.env.JEST_WORKER_ID,
        isIntegrationTest: this._getBoolean('INTEGRATION_TEST', false),
        enableDebug: this._getBoolean('DEBUG', false),
        enableScreenshots: this._getBoolean('ENABLE_ERROR_SCREENSHOTS', false)
      }
    };

    return this._config;
  }

  /**
   * Type-safe configuration getters
   * These methods encapsulate parsing and validation logic
   */
  _getString(key, defaultValue = undefined) {
    const value = process.env[key];
    if (value === undefined) return defaultValue;
    return value.trim();
  }

  _getNumber(key, defaultValue = undefined) {
    const value = process.env[key];
    if (value === undefined) return defaultValue;
    const parsed = parseInt(value, 10);
    if (isNaN(parsed)) {
      console.warn(`⚠️ Invalid number for ${key}: ${value}, using default: ${defaultValue}`);
      return defaultValue;
    }
    return parsed;
  }

  _getBoolean(key, defaultValue = false) {
    const value = process.env[key];
    if (value === undefined) return defaultValue;
    return value.toLowerCase() === 'true';
  }

  /**
   * Public Interface - Simple methods that hide complexity
   */

  get(path = null) {
    const config = this._buildConfig();
    if (!path) return config;
    
    return path.split('.').reduce((obj, key) => obj?.[key], config);
  }

  isProduction() {
    return this.get('nodeEnv') === 'production';
  }

  isDevelopment() {
    return this.get('nodeEnv') === 'development';
  }

  isTest() {
    return this.get('development.isTest');
  }

  isServerless() {
    const cloud = this.get('cloud');
    return cloud.isCloudFunction || cloud.isAppEngine || cloud.isNetlify;
  }

  getPaymentConfig() {
    const payments = this.get('payments');
    return {
      enabled: payments.enabled,
      devMode: payments.devMode,
      required: payments.required,
      publishableKey: payments.stripePublishableKey
    };
  }

  getDatabaseConfig() {
    return this.get('database');
  }

  /**
   * Validation with clear error messages
   * Deep module: complex validation logic with simple interface
   */
  validate() {
    const config = this._buildConfig();
    const errors = [];

    // Production validations
    if (this.isProduction()) {
      if (!config.openai.apiKey) {
        errors.push('OPENAI_API_KEY is required in production');
      }
      
      if (config.database.enabled && !config.database.password) {
        errors.push('DB_PASSWORD is required when database is enabled in production');
      }
    }

    // Port validations
    if (config.port < 1 || config.port > 65535) {
      errors.push(`Invalid PORT: ${config.port}. Must be between 1 and 65535`);
    }

    if (errors.length > 0) {
      const message = 'Configuration validation failed:\n' + errors.map(e => `  - ${e}`).join('\n');
      throw new Error(message);
    }

    return true;
  }

  /**
   * Debug information for troubleshooting
   */
  getDebugInfo() {
    const config = this._buildConfig();
    return {
      nodeEnv: config.nodeEnv,
      port: config.port,
      hasOpenAIKey: !!config.openai.apiKey,
      databaseEnabled: config.database.enabled,
      paymentsEnabled: config.payments.enabled,
      isServerless: this.isServerless(),
      initialized: this._initialized
    };
  }
}

// Export singleton instance
const configuration = new Configuration();

module.exports = configuration;


