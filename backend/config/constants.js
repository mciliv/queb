/**
 * APPLICATION CONFIGURATION CONSTANTS
 * Purpose: Centralized configuration for all application components
 * Benefits: Single source of truth for constants, easier maintenance, clearer intentions
 */

// ==================== SERVER CONFIGURATION ====================
const SERVER_CONFIG = {
  // Port configuration
  DEFAULT_HTTP_PORT: 8080,
  DEFAULT_HTTPS_PORT: 3001,
  LIVERELOAD_PORT: 35730,
  
  // Timeout settings
  REQUEST_TIMEOUT: 45000,  // 45 seconds for API requests
  AI_TIMEOUT: 30000,       // 30 seconds for AI API calls
  DATABASE_TIMEOUT: 2000,  // 2 seconds for database connections
  
  // Development settings
  MOBILE_IP: '172.20.10.4'  // Default mobile development IP
};

// ==================== AI SERVICE CONFIGURATION ====================
const AI_CONFIG = {
  // OpenAI settings
  MODEL: 'gpt-4o',
  MAX_TOKENS: 1000,
  TEMPERATURE: 0.1,
  MAX_RETRIES: 2,
  
  // Image processing
  IMAGE_DETAIL: 'high',
  CROP_SIZE_DEFAULT: 200
};

// ==================== DATABASE CONFIGURATION ====================
const DATABASE_CONFIG = {
  // Connection defaults
  DEFAULT_HOST: 'localhost',
  DEFAULT_PORT: 5432,
  DEFAULT_DATABASE: 'mol_users',
  DEFAULT_USER: 'mol_user',
  DEFAULT_PASSWORD: 'mol_password',
  
  // Pool settings
  MAX_CONNECTIONS: 20,
  IDLE_TIMEOUT: 30000,
  CONNECTION_TIMEOUT: 2000
};

// ==================== FILE SYSTEM CONFIGURATION ====================
const FILE_CONFIG = {
  // Directories
  SDF_DIRECTORY: 'data/sdf_files',
  CHEMISTRY_PATH: 'chemistry/processors',
  
  // File extensions
  WATCHED_EXTENSIONS: ['.js', '.html', '.css'],
  
  // Processing limits
  MAX_FILE_SIZE: 10485760,  // 10MB in bytes
  MAX_ERRORS_STORED: 50
};

// ==================== UI CONFIGURATION ====================
const UI_CONFIG = {
  // Molecular visualization
  SPHERE_SCALE: 0.8,  // van der Waals radii scale for 3D molecules
  BACKGROUND_COLOR: 'black',
  
  // Animation and timing
  RENDER_DELAY: 100,    // Stagger molecular rendering (ms)
  RESIZE_DELAY: 50,     // Viewer resize delay (ms)
  ERROR_DISPLAY_TIME: 5000,  // Auto-remove errors after 5 seconds
  
  // Layout
  PROCESSING_Z_INDEX: 1000
};

// ==================== ERROR HANDLING CONFIGURATION ====================
const ERROR_CONFIG = {
  // HTTP Status codes for different error types
  STATUS_CODES: {
    VALIDATION_ERROR: 400,
    UNAUTHORIZED: 401,
    TIMEOUT: 408,
    RATE_LIMITED: 429,
    INTERNAL_ERROR: 500,
    SERVICE_UNAVAILABLE: 503
  },
  
  // Error prefixes for logging
  PREFIXES: {
    ERROR: '🚨',
    SUCCESS: '✅',
    WARNING: '⚠️',
    INFO: '💡'
  }
};

// ==================== CHEMISTRY CONFIGURATION ====================
const CHEMISTRY_CONFIG = {
  // Validation patterns
  MOLECULAR_FORMULA_PATTERN: /^[A-Z][0-9]*([A-Z][0-9]*)*$/,
  
  // Processing limits
  MAX_SMILES_LENGTH: 1000,
  SMILES_BATCH_SIZE: 10,
  
  // Invalid SMILES indicators
  INVALID_SMILES: ['N/A', '', null, undefined],
  
  // Chemical representation formats
  SUPPORTED_FORMATS: ['SMILES', 'SDF', 'MOL']
};

// ==================== ENVIRONMENT DETECTION ====================
const ENVIRONMENT = {
  isDevelopment: () => process.env.NODE_ENV === 'development',
  isProduction: () => process.env.NODE_ENV === 'production',
  isLocalhost: () => ['localhost', '127.0.0.1'].includes(location?.hostname),
  
  // Feature flags based on environment
  FEATURES: {
    AUTO_ENABLE_DEV_MODE: true,  // Auto-enable dev mode on localhost
    VERBOSE_LOGGING: process.env.NODE_ENV === 'development',
    ERROR_STACK_TRACES: process.env.NODE_ENV === 'development'
  }
};

// ==================== EXPORTS ====================
module.exports = {
  SERVER_CONFIG,
  AI_CONFIG,
  DATABASE_CONFIG,
  FILE_CONFIG,
  UI_CONFIG,
  ERROR_CONFIG,
  CHEMISTRY_CONFIG,
  ENVIRONMENT
};

// For frontend usage (if imported via script tag)
if (typeof window !== 'undefined') {
  window.APP_CONFIG = {
    SERVER_CONFIG,
    AI_CONFIG,
    UI_CONFIG,
    ERROR_CONFIG,
    CHEMISTRY_CONFIG,
    ENVIRONMENT
  };
}