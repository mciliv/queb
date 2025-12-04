// test/suites/integration/fixtures/index.js - Main fixtures export
const fixtures = require('./fixtures');
const utils = require('./utils');

// Export all fixtures and utilities
module.exports = {
  ...fixtures,
  ...utils,
};

