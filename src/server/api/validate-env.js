// Environment validation for local/dev startup
// Validates required environment variables and configuration

function validateLocalDevEnv(config, logger) {
  // Only validate in local/dev environments, not in production or test
  if (config.isProduction() || config.isTest()) {
    return; // Skip validation in production/test
  }

  const errors = [];
  const warnings = [];

  // Check OpenAI configuration
  const openaiApiKey = config.get('openai.apiKey');
  const openaiModel = config.get('openai.model');

  if (!openaiApiKey) {
    errors.push('OPENAI_API_KEY is required for local development');
  }

  if (!openaiModel) {
    errors.push('OPENAI_MODEL is required for local development');
  } else if (openaiModel && !['gpt-4o', 'gpt-4', 'gpt-3.5-turbo'].includes(openaiModel)) {
    warnings.push(`OPENAI_MODEL=${openaiModel} may not be supported. Recommended: gpt-4o, gpt-4, or gpt-3.5-turbo`);
  }

  // Log warnings (non-fatal)
  if (warnings.length > 0) {
    warnings.forEach(warning => logger.warn(`⚠️  ${warning}`));
  }

  // Throw errors (fatal)
  if (errors.length > 0) {
    const errorMessage = 'Environment validation failed:\n' + errors.map(e => `  - ${e}`).join('\n');
    logger.error('❌ ' + errorMessage);
    throw new Error(errorMessage);
  }

  // Silent on success - only report errors/warnings
}

module.exports = {
  validateLocalDevEnv
};
