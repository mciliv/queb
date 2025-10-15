const path = require('path');
const tomlConfig = require('../../config');

class Configuration {
  constructor() {
    this._config = null;
    this._initialized = false;
  }

  _buildConfigWithTypeSafetyAndValidation() {
    if (this._config) return this._config;
    if (!(tomlConfig && typeof tomlConfig.getConfig === 'function')) {
      throw new Error('Configuration error: TOML-based configuration is required but not found. Ensure config/config.toml exists.');
    }

    const cfg = tomlConfig.getConfig();
    this._config = {
<<<<<<< Updated upstream
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
=======
      nodeEnv: cfg.environment?.nodeEnv ?? 'development',
      port: cfg.server?.port ?? 8080,
      openai: cfg.openai ?? {},
      database: cfg.database ?? {},
      payments: cfg.payments ?? {},
      cloud: cfg.cloud ?? {},
      ssl: cfg.ssl ?? {},
      development: Object.assign({
        isTest: (cfg.environment?.nodeEnv === 'test')
      }, cfg.development || {}),
      enableAllFeatures: cfg.enableAllFeatures === true || (cfg.development && cfg.development.enableAllFeatures === true) || (cfg.features && cfg.features.enableAll === true)
>>>>>>> Stashed changes
    };

    // Consolidation logic: If enableAllFeatures is true, override individual service enables
    const enableAllFeatures = this._config.enableAllFeatures === true;
    if (enableAllFeatures) {
      if (this._config.database) this._config.database.enabled = true;
      if (this._config.payments) this._config.payments.enabled = true;
    }
    // Mark as initialized once config is built
    this._initialized = true;
    return this._config;
  }

  // Env-parsing helpers removed: TOML is the single source of truth

  get(path = null) {
    const config = this.getConfig();
    if (!path) return config;
    
    return path.split('.').reduce((obj, key) => obj?.[key], config);
  }

  // Public single source of truth for configuration
  getConfig() {
    return this._buildConfigWithTypeSafetyAndValidation();
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

  validate() {
    const config = this._buildConfigWithTypeSafetyAndValidation();
    const validationErrors = [];

    const isProductionEnvironment = this.isProduction();
    if (isProductionEnvironment) {
      const openAIKeyMissing = !config.openai.apiKey;
      if (openAIKeyMissing) {
        validationErrors.push('OPENAI_API_KEY is required in production');
      }
      
      const databasePasswordRequiredButMissing = config.database.enabled && !config.database.password;
      if (databasePasswordRequiredButMissing) {
        validationErrors.push('DB_PASSWORD is required when database is enabled in production');
      }
    }

    const portIsInvalid = config.port < 1 || config.port > 65535;
    if (portIsInvalid) {
      validationErrors.push(`Invalid PORT: ${config.port}. Must be between 1 and 65535`);
    }

    const hasValidationErrors = validationErrors.length > 0;
    if (hasValidationErrors) {
      const errorMessage = 'Configuration validation failed:\n' + validationErrors.map(e => `  - ${e}`).join('\n');
      throw new Error(errorMessage);
    }

    return true;
  }

  getDebugInfo() {
    const config = this._buildConfigWithTypeSafetyAndValidation();
    return {
      nodeEnv: config.nodeEnv,
      port: config.port,
      hasOpenAIKey: !!config.openai.apiKey,
      databaseEnabled: config.database.enabled,
      paymentsEnabled: config.payments.enabled,
      allFeaturesEnabled: !!config.enableAllFeatures,
      isServerless: this.isServerless(),
      initialized: this._initialized
    };
  }
}

const configurationSingleton = new Configuration();

module.exports = configurationSingleton;
