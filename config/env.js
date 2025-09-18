// Backend Environment Configuration - now consolidated in project.js
// This file is deprecated, use project.js configuration instead

const project = require('./project');

module.exports = {
  // Environment
  NODE_ENV: project.env.NODE_ENV,
  PORT: project.env.PORT,

  // API
  OPENAI_API_KEY: project.api.OPENAI_API_KEY,

  // Payments
  PAYMENTS_ENABLED: project.payments.enabled,
  PAYMENTS_DEV_MODE: project.payments.devMode,
  PAYMENTS_REQUIRED: project.payments.required,

  // Database
  DB_ENABLED: project.database.enabled,
  DB_HOST: project.database.host,
  DB_PORT: project.database.port,
  DB_NAME: project.database.name,
  DB_USER: project.database.user,
  DB_PASSWORD: project.database.password,

  // Cloud
  GOOGLE_CLOUD_PROJECT: project.cloud.googleCloudProject,
  FUNCTION_NAME: project.cloud.functionName,
  FUNCTION_TARGET: project.cloud.functionTarget,
  K_SERVICE: project.cloud.kService,

  // Test
  INTEGRATION_TEST: project.test.integration,

  // SSL
  SSL_CERT_PATH: project.ssl.certPath,
  SSL_KEY_PATH: project.ssl.keyPath,

  // Helper functions
  getPaymentConfig: () => project.helpers.getPaymentConfig(),
  validateConfig: () => project.validate()
}; 