// Modern DotenvX - Next-gen environment management
// Features: encryption, multi-environment, secure by default

const { config } = require('@dotenvx/dotenvx');

// Load with encryption and multi-environment support
function loadDotenvX() {
  // Load configuration with encryption support
  const result = config({
    // Load multiple files in priority order
    path: [
      '.env.defaults',     // Default values
      '.env',              // Project configuration
      '.env.local',        // Local overrides
      '.env.production',   // Production overrides (if NODE_ENV=production)
    ],

    // Enable encryption for sensitive values
    encryption: true,

    // Environment-specific loading
    env: process.env.NODE_ENV || 'development',

    // Override existing environment variables
    override: true,

    // Debug mode for troubleshooting
    debug: process.env.DEBUG === 'true'
  });

  if (result.error) {
    console.error('❌ DotenvX configuration error:', result.error.message);
    process.exit(1);
  }

  console.log('✅ Loaded configuration with DotenvX');
  if (result.debug) {
    console.log('📊 Loaded variables:', Object.keys(result.parsed || {}).length);
  }

  return result.parsed;
}

// Usage examples:
//
// 1. Basic usage:
//    node -r ./config/dotenvx-config.js your-app.js
//
// 2. With encryption:
//    npx dotenvx encrypt .env.local
//    npx dotenvx run -- node your-app.js
//
// 3. Multi-environment:
//    NODE_ENV=production node -r ./config/dotenvx-config.js your-app.js
//
// 4. Debug mode:
//    DEBUG=true node -r ./config/dotenvx-config.js your-app.js
//
// 5. CI/CD with secrets:
//    npx dotenvx run -- node your-app.js
//
// 6. Compare environments:
//    npx dotenvx compare .env .env.production
//
// 7. Generate .env.example:
//    npx dotenvx get > .env.example

module.exports = { loadDotenvX };
