// Client-side application configuration
// No server-side dependencies for browser compatibility

export const PAYMENT_CONFIG = {
  enabled: process.env.NODE_ENV === 'development',
  devMode: process.env.NODE_ENV === 'development',
  required: false,
  effectiveEnabled: process.env.NODE_ENV === 'development'
};

export const APP_CONFIG = {
  // Client-side configuration
  environment: process.env.NODE_ENV || 'development',
  version: '1.0.0'
};
