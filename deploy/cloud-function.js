// Cloud Functions 2nd Gen Entry Point
// This is a simplified entry point specifically for Cloud Functions deployment

const express = require('express');
const cors = require('cors');
const path = require('path');

// Load configuration
const config = require('./config/env');

// Create Express app
const app = express();

// Middleware
app.use(cors());
app.use(express.json({ limit: '10mb' }));
app.use(express.urlencoded({ extended: true, limit: '10mb' }));

// Serve static files
app.use(express.static(path.join(__dirname, 'src/client')));
app.use(express.static(path.join(__dirname, 'public')));

// Health check endpoint
app.get('/health', (req, res) => {
  res.json({ status: 'ok', timestamp: new Date().toISOString() });
});

// Main API endpoint
app.post('/molecularAnalysis', async (req, res) => {
  try {
    // Import the main server logic
    const mainApp = require('./src/server/api/server');
    
    // Handle the request
    if (typeof mainApp === 'function') {
      return mainApp(req, res);
    } else if (mainApp.molecularAnalysis) {
      return mainApp.molecularAnalysis(req, res);
    } else {
      res.status(500).json({ error: 'Function not available' });
    }
  } catch (error) {
    console.error('Cloud Function error:', error);
    res.status(500).json({ error: error.message });
  }
});

// Catch-all handler for other routes
app.all('*', (req, res) => {
  res.status(404).json({ error: 'Not found' });
});

// Start server on PORT environment variable
const port = process.env.PORT || 8080;
const server = app.listen(port, '0.0.0.0', () => {
  console.log(`✅ Cloud Function server started on port ${port}`);
});

// Handle server errors
server.on('error', (error) => {
  console.error('❌ Server error:', error);
  process.exit(1);
});

// Export for Cloud Functions
exports.molecularAnalysis = app;

