// Environment Configuration
// Loads environment variables from .env file
// Note: .env file is hidden from assistant for security

const path = require('path');
require('dotenv').config({ 
  path: path.resolve(process.cwd(), '.env'),
  debug: false, // Set to true for debugging env loading
  quiet: true  // Suppress all dotenv logs
});

const config = {
  // Node Environment
  NODE_ENV: process.env.NODE_ENV || 'development',
  
  // Server Configuration
  PORT: parseInt(process.env.PORT) || 8080,
  
  // OpenAI API Configuration
  OPENAI_API_KEY: process.env.OPENAI_API_KEY,
  
  // Payment Configuration
  PAYMENTS_ENABLED: process.env.PAYMENTS_ENABLED === 'false',
  PAYMENTS_DEV_MODE: process.env.PAYMENTS_DEV_MODE === 'true',
  PAYMENTS_REQUIRED: process.env.PAYMENTS_REQUIRED === 'false',
  
  // Database Configuration (PostgreSQL)
  DB_HOST: process.env.DB_HOST || 'localhost',
  DB_PORT: parseInt(process.env.DB_PORT) || 5432,
  DB_NAME: process.env.DB_NAME || 'mol_users',
  DB_USER: process.env.DB_USER || 'mol_user',
  DB_PASSWORD: process.env.DB_PASSWORD || 'mol_password',
  
  // Cloud Configuration (optional)
  GOOGLE_CLOUD_PROJECT: process.env.GOOGLE_CLOUD_PROJECT,
  FUNCTION_NAME: process.env.FUNCTION_NAME,
  FUNCTION_TARGET: process.env.FUNCTION_TARGET,
  K_SERVICE: process.env.K_SERVICE,
  
  // Netlify Configuration (optional)
  NETLIFY: process.env.NETLIFY,
  
  // Test Configuration
  INTEGRATION_TEST: process.env.INTEGRATION_TEST === 'true',
  
  // SSL/HTTPS Configuration (optional)
  SSL_CERT_PATH: process.env.SSL_CERT_PATH,
  SSL_KEY_PATH: process.env.SSL_KEY_PATH,
};

// Payment configuration helper
const getPaymentConfig = () => {
  const isDev = config.NODE_ENV === 'development';
  const isDevDomain = process.env.HOSTNAME === 'localhost' || process.env.HOSTNAME === '127.0.0.1';
  
  return {
    enabled: config.PAYMENTS_ENABLED || false,
    devMode: config.PAYMENTS_DEV_MODE || isDev || isDevDomain,
    required: config.PAYMENTS_REQUIRED || false,
    // Auto-disable in development unless explicitly enabled
    effectiveEnabled: config.PAYMENTS_ENABLED === true || (isDev && config.PAYMENTS_DEV_MODE !== false)
  };
};

// Validation
const validateConfig = () => {
  const errors = [];
  
  // Required for production
  if (config.NODE_ENV === 'production' && !config.OPENAI_API_KEY) {
    errors.push('OPENAI_API_KEY is required in production');
  }
  
  if (errors.length > 0) {
    throw new Error(`Configuration errors:\n${errors.join('\n')}`);
  }
};

// Export validated config
module.exports = {
  ...config,
  getPaymentConfig,
  validateConfig,
}; 