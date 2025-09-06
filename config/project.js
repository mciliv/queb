// Consolidated Project Configuration
// Contains all essential project settings and configuration

const path = require('path');

// Load environment variables
require('dotenv').config({
  path: path.resolve(process.cwd(), '.env'),
  debug: false,
  quiet: true
});

module.exports = {
  // Basic Project Information
  name: 'molecular-analysis',
  description: 'Molecular analysis app with AI-powered chemical identification',
  version: process.env.npm_package_version || '1.0.0',

  // Environment Configuration
  env: {
    NODE_ENV: process.env.NODE_ENV || 'development',
    PORT: parseInt(process.env.PORT) || 8080,
    HOSTNAME: process.env.HOSTNAME || 'localhost'
  },

  // API Configuration
  api: {
    OPENAI_API_KEY: process.env.OPENAI_API_KEY,
    PYTHON_VERSION: process.env.PYTHON_VERSION || '3.12.9'
  },

  // Database Configuration
  database: {
    enabled: process.env.DB_ENABLED !== 'false',
    host: process.env.DB_HOST || 'localhost',
    port: parseInt(process.env.DB_PORT) || 5432,
    name: process.env.DB_NAME || 'mol_users',
    user: process.env.DB_USER || 'mol_user',
    password: process.env.DB_PASSWORD || 'mol_password'
  },

  // Payment Configuration
  payments: {
    enabled: process.env.PAYMENTS_ENABLED === 'true',
    devMode: process.env.PAYMENTS_DEV_MODE === 'true',
    required: process.env.PAYMENTS_REQUIRED === 'true'
  },

  // Cloud Configuration
  cloud: {
    domain: process.env.DOMAIN_NAME || '',
    dnsZone: process.env.DNS_ZONE_NAME || '',
    region: process.env.REGION || 'us-central1',
    functionName: process.env.FUNCTION_NAME || 'molecular-analysis',
    googleCloudProject: process.env.GOOGLE_CLOUD_PROJECT,
    functionTarget: process.env.FUNCTION_TARGET,
    kService: process.env.K_SERVICE
  },

  // SSL Configuration
  ssl: {
    certPath: process.env.SSL_CERT_PATH,
    keyPath: process.env.SSL_KEY_PATH
  },

  // Test Configuration
  test: {
    integration: process.env.INTEGRATION_TEST === 'true',
    visual: process.env.REACT_APP_RUN_VISUAL_TESTS === 'true'
  },

  // Version Configuration
  version: {
    defaultVersion: 'react', // 'vanilla' | 'react'
    allowToggle: false,
    showToggleButton: false,
    persistChoice: false,
    performanceMode: false
  },

  // Frontend Constants
  constants: {
    testMolecules: {
      water: { name: 'Water', smiles: 'O' },
      ethanol: { name: 'Ethanol', smiles: 'CCO' },
      caffeine: { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' }
    },
    visualTests: [
      {
        label: 'Coffee',
        smilesList: ['CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'O', 'CCO']
      },
      {
        label: 'Red wine',
        smilesList: ['O', 'CCO', 'CC(=O)O', 'C1=CC=CC=C1']
      },
      {
        label: 'Apple',
        smilesList: ['O', 'CC(=O)O']
      }
    ]
  },

  // Helper Functions
  helpers: {
    getPaymentConfig() {
      const isDev = this.env.NODE_ENV === 'development';
      const isDevDomain = this.env.HOSTNAME === 'localhost' || this.env.HOSTNAME === '127.0.0.1';

      return {
        enabled: this.payments.enabled,
        devMode: this.payments.devMode || isDev || isDevDomain,
        required: this.payments.required,
        effectiveEnabled: this.payments.enabled || (isDev && this.payments.devMode !== false)
      };
    },

    getVersionConfig() {
      // Environment-specific overrides
      const envConfig = this.env.NODE_ENV === 'production'
        ? { allowToggle: false, showToggleButton: false, performanceMode: false }
        : { allowToggle: false, showToggleButton: false, performanceMode: false };

      return {
        ...this.version,
        ...envConfig
      };
    },

    getSmilesNameMap() {
      return Object.values(this.constants.testMolecules).reduce((acc, m) => {
        acc[m.smiles] = m.name;
        return acc;
      }, {});
    }
  },

  // Validation
  validate() {
    const errors = [];

    if (this.env.NODE_ENV === 'production' && !this.api.OPENAI_API_KEY) {
      errors.push('OPENAI_API_KEY is required in production');
    }

    if (errors.length > 0) {
      throw new Error(`Configuration errors:\n${errors.join('\n')}`);
    }

    return true;
  }
};

// Auto-validate on load
module.exports.validate();
