// Cloud Function entry point
const app = require('./backend/api/server');

// Export for Google Cloud Functions
exports.molecularAnalysis = app;
