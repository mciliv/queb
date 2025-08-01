// Environment Loader - Assumes all environment variables are already loaded
// This file doesn't depend on .env directly, works with any environment source

const envLoader = {
  // Server Configuration
  PORT: process.env.PORT || 8080,
  NODE_ENV: process.env.NODE_ENV || 'development',
  
  // API Keys (assumes already loaded from .env or other source)
  OPENAI_API_KEY: process.env.OPENAI_API_KEY,
  
  // Database Configuration
  DB_HOST: process.env.DB_HOST || 'localhost',
  DB_PORT: parseInt(process.env.DB_PORT) || 5432,
  DB_NAME: process.env.DB_NAME || 'mol_users',
  DB_USER: process.env.DB_USER || 'mol_user',
  DB_PASSWORD: process.env.DB_PASSWORD || 'mol_password',
  
  // Payment Configuration
  STRIPE_PUBLISHABLE_KEY: process.env.STRIPE_PUBLISHABLE_KEY,
  STRIPE_SECRET_KEY: process.env.STRIPE_SECRET_KEY,
  
  // Cloud Configuration
  GOOGLE_CLOUD_PROJECT: process.env.GOOGLE_CLOUD_PROJECT,
  FUNCTION_NAME: process.env.FUNCTION_NAME,
  
  // SSL Configuration
  SSL_CERT_PATH: process.env.SSL_CERT_PATH,
  SSL_KEY_PATH: process.env.SSL_KEY_PATH,
  
  // Get all environment variables
  getAll() {
    return { ...this };
  },
  
  // Check if required variables are present
  validate() {
    const missing = [];
    
    if (!this.OPENAI_API_KEY) {
      missing.push('OPENAI_API_KEY');
    }
    
    return {
      isValid: missing.length === 0,
      missing: missing
    };
  },
  
  // Log configuration status (without sensitive values)
  logStatus() {
    const validation = this.validate();
    console.log('🔧 Environment Configuration:');
    console.log(`   NODE_ENV: ${this.NODE_ENV}`);
    console.log(`   PORT: ${this.PORT}`);
    console.log(`   DB_HOST: ${this.DB_HOST}`);
    console.log(`   OPENAI_API_KEY: ${this.OPENAI_API_KEY ? 'Present' : 'Missing'}`);
    
    if (!validation.isValid) {
      console.log(`⚠️  Missing variables: ${validation.missing.join(', ')}`);
    } else {
      console.log('✅ All required variables present');
    }
  }
};

module.exports = envLoader;