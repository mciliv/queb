// Toggle one option below by uncommenting it. Keep others commented.

// 1) Cloud Functions (default)
const app = require('./backend/api/server');
exports.molecularAnalysis = app;

// 2) Cursor IDE preview (local server)
// const app = require('./backend/api/server');
// server.js auto-starts local HTTP/HTTPS when not serverless
// module.exports = app; // optional
