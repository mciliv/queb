// Application configuration - now consolidated in project.js
// This file is deprecated, use project.js configuration instead

import project from '../../config/project.js';

export const PAYMENT_CONFIG = project.helpers.getPaymentConfig();
export const APP_CONFIG = {
  // All configuration moved to project.js
};
