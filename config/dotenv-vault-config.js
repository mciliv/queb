// Dotenv Vault - Enterprise-grade secrets management
// Perfect for team environments and production deployment

const { config } = require('dotenv');
const { config: vaultConfig } = require('dotenv-vault-core');

// Load environment with vault support
function loadVaultConfig() {
  // Try vault first (for production/team environments)
  try {
    vaultConfig();
    console.log('✅ Loaded configuration from Dotenv Vault');
  } catch (vaultError) {
    // Fallback to local .env files
    console.log('⚠️  Vault not available, using local .env files');

    // Load in hierarchical order
    config({ path: '.env.defaults' }); // Defaults
    config({ path: '.env' });          // Project config
    config({ path: '.env.local' });    // Local overrides

    console.log('✅ Loaded configuration from local files');
  }
}

// Usage examples:
//
// 1. Development (local files):
//    node your-app.js
//
// 2. Production (vault):
//    DOTENV_KEY=dotenv://:key_1234567890abcdef@dotenv.org/vault/.env.vault?environment=production node your-app.js
//
// 3. Team collaboration:
//    - Run: npx dotenv-vault push
//    - Team members: npx dotenv-vault pull
//
// 4. CI/CD:
//    - Store DOTENV_KEY as environment variable
//    - Vault handles encryption/decryption automatically

module.exports = { loadVaultConfig };
