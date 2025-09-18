// Configuration loader utility
import { PAYMENT_CONFIG, APP_CONFIG } from '../config/app-config.js';
import { VALIDATION_PATTERNS } from '../config/validation-patterns.js';

/**
 * Centralized configuration loader
 * Provides a single point of access for all application configurations
 */
export const ConfigLoader = {
  // Payment configuration
  getPaymentConfig: () => PAYMENT_CONFIG,
  
  // App configuration
  getAppConfig: () => APP_CONFIG,
  
  // Validation patterns
  getValidationPatterns: () => VALIDATION_PATTERNS,
  
  // Combined getter for convenience
  getAllConfigs: () => ({
    payment: PAYMENT_CONFIG,
    app: APP_CONFIG,
    validation: VALIDATION_PATTERNS
  })
};

// Direct exports for backward compatibility
export { PAYMENT_CONFIG, APP_CONFIG, VALIDATION_PATTERNS };
