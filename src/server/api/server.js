const express = require('express');
const cors = require('cors');
const path = require('path');
const fs = require('fs');
const createError = require('http-errors');
const { createContainer } = require('../../core/services');

/**
 * Create and configure Express application with injected services
 * @param {ServiceContainer} container - Configured DI container
 * @returns {express.Application} Configured Express app
 */
async function createApp(container) {
  const app = express();

  // Get core services
  const config = await container.get('config');
  const logger = await container.get('logger');
  const errorHandler = await container.get('errorHandler');

  // Middleware
  app.use(cors());
  app.use(express.json({ limit: '50mb' }));
  app.use(express.urlencoded({ extended: true, limit: '50mb' }));

  // Define serveIndexHtml function for dynamic HTML serving
  const distPath = path.join(__dirname, '..', '..', 'client', 'dist');
  const serveIndexHtml = (req, res) => {
    try {
      const indexPath = path.join(distPath, 'index.html');
      let htmlContent = fs.readFileSync(indexPath, 'utf8');

      // Dynamic token replacement for SEO
      const baseUrl = `${req.protocol}://${req.get('host')}`;
      const cacheBust = Date.now();

      htmlContent = htmlContent
        .replace(/\{\{TITLE\}\}/g, 'Queb - Molecular Analysis')
        .replace(/\{\{DESCRIPTION\}\}/g, '3D visualization of chemical contents by camera or text input.')
        .replace(/\{\{CANONICAL\}\}/g, `${baseUrl}${req.path}`)
        .replace(/\{\{OG_IMAGE\}\}/g, '/images/favicon.svg')
        .replace(/\{\{CACHE_BUST\}\}/g, cacheBust);

      res.type('html').send(htmlContent);
    } catch (error) {
      logger.error('Failed to serve index.html:', error);
      res.status(500).send('Server error');
    }
  };

  // Root route serves the frontend with token replacement (must be BEFORE static middleware)
  app.get('/', serveIndexHtml);

  // Serve static files from public directory (images, manifest, etc)
  const publicPath = path.join(__dirname, '..', '..', '..', 'public');
  app.use(express.static(publicPath));

  // Serve /assets from client assets directory (CSS, icons)
  const assetsPath = path.join(__dirname, '..', '..', 'client', 'assets');
  app.use('/assets', express.static(assetsPath));

  // Serve /dist from client dist directory (bundle.js, bundle.css, etc)
  app.use('/dist', express.static(distPath));

  // Serve other dist files at root (bundle.js, bundle.css without /dist prefix)
  // index: false prevents this from serving index.html
  app.use(express.static(distPath, { index: false }));

  // Request logging middleware
  app.use((req, res, next) => {
    const start = Date.now();
    const originalSend = res.send;
    
    res.send = function(data) {
      res.send = originalSend;
      const duration = Date.now() - start;
      logger.info(`${req.method} ${req.path} - ${res.statusCode} (${duration}ms)`);
      return res.send(data);
    };
    
    next();
  });
  
  // Health check
  app.get('/api/health', (req, res) => {
    res.json({ 
      status: 'ok', 
      environment: config.get('nodeEnv'),
      version: process.env.npm_package_version || '1.0.0'
    });
  });

  // Frontend error logging endpoint (dev only)
  app.post('/api/log-error', (req, res) => {
    const { type, message, timestamp, source, stack, location, url } = req.body;
    
    // Log to server logger
    const logMessage = `[${source || 'frontend'}] ${type || 'error'}: ${message || 'Unknown error'}`;
    const logData = {
      type,
      message,
      timestamp,
      source,
      location,
      url,
      stack: stack ? stack.substring(0, 500) : undefined // Truncate long stacks
    };

    if (type === 'error') {
      logger.error(logMessage, logData);
    } else if (type === 'warn') {
      logger.warn(logMessage, logData);
    } else {
      logger.info(logMessage, logData);
    }

    // Always return success to avoid cluttering console with failed requests
    res.status(200).json({ received: true });
  });
  
  await setupChemicalPredictionRoutes(app, container);

  // Setup user routes if database is enabled
  if (config.get('database.enabled')) {
    const { setupUserRoutes } = require('./user-routes');
    await setupUserRoutes(app, container);
  }

  // SPA catch-all - serve index.html for non-API, non-file requests
  app.get('*', (req, res, next) => {
    // Skip if it's an API route or file request
    if (req.path.startsWith('/api/') ||
        req.path.startsWith('/sdf_files/') ||
        req.path.includes('.')) {
      return next();
    }

    serveIndexHtml(req, res); // For SPA client-side routing
  });

  // 404 handler - catch unmatched routes and pass to error handler
  app.use((req, res, next) => {
    // For API routes, create a 404 error with JSON-friendly message
    if (req.path.startsWith('/api/')) {
      return next(createError(404, 'API endpoint not found'));
    }
    // For other routes, create a generic 404
    next(createError(404, 'Not found'));
  });

  // Error handling middleware
  app.use(async (err, req, res, next) => {
    // Respect HTTP status code from http-errors (defaults to 500)
    const status = err.status || err.statusCode || 500;
    
    // If error already has structured properties (from ErrorHandler), use them
    // Otherwise, handle the raw error
    let handled;
    if (err.recoverable !== undefined && err.recovery !== undefined) {
      handled = {
        message: err.message,
        code: err.code,
        recoverable: err.recoverable,
        recovery: err.recovery,
        timestamp: err.timestamp,
        context: err.context
      };
    } else {
      handled = await errorHandler.handle(err, {
        path: req.path,
        method: req.method,
        ip: req.ip
      });
    }

    // Only log errors (not 404s) to avoid noise
    if (status >= 500) {
      logger.error(`API Error: ${handled.message}`, {
        error: handled,
        request: {
          path: req.path,
          method: req.method,
          body: req.body
        }
      });
    }

    const isDev = process.env.NODE_ENV === 'dev';

    // For API routes, always return JSON
    if (req.path.startsWith('/api/')) {
      return res.status(status).json({
        error: handled.message || err.message || 'Internal server error',
        code: handled.code,
        recoverable: handled.recoverable,
        recovery: handled.recovery,
        timestamp: handled.timestamp,
        ...(isDev && {
          stack: err.stack,
          context: handled.context,
          originalMessage: err.message
        })
      });
    }

    // For non-API routes, return appropriate response based on status
    if (status === 404) {
      return res.status(404).send('Not found');
    }
    
    res.status(status).json({
      error: handled.message || err.message || 'Internal server error',
      ...(isDev && {
        stack: err.stack,
        context: handled.context
      })
    });
  });
  
  return app;
}

async function setupChemicalPredictionRoutes(app, container) {
  const structuralizer = await container.get('structuralizer');
  const molecularProcessor = await container.get('molecularProcessor');
  const logger = await container.get('logger');

  app.post('/api/structuralize', async (req, res, next) => {
    try {
      const { text, lookupMode = 'GPT-5' } = req.body;

      if (!text || typeof text !== 'string') {
        return res.status(400).json({
          error: 'Missing or invalid "text" parameter'
        });
      }

      if (typeof lookupMode !== 'string') {
        return res.status(400).json({
          error: 'Invalid "lookupMode" parameter'
        });
      }

      logger.info(`Text prediction request: "${text}" (mode: ${lookupMode})`);

      const result = await structuralizer.chemicals({
        object: text.trim(),
        lookupMode
      });

      res.json(result);
    } catch (error) {
      next(error);
    }
  });
  
  app.post('/api/structuralize-image', async (req, res, next) => {
    try {
      const { imageBase64, x, y } = req.body;
      
      if (!imageBase64) {
        return res.status(400).json({ 
          error: 'Missing imageBase64 parameter' 
        });
      }
      
      logger.info('Image prediction request', {
        hasCoordinates: x !== undefined && y !== undefined
      });
      
      const result = await structuralizer.chemicals({
        imageBase64,
        x,
        y
      });
      
      res.json(result);
    } catch (error) {
      next(error);
    }
  });
  
  // SDF generation endpoint
  app.post('/api/generate-sdfs', async (req, res, next) => {
    try {
      const { smiles, overwrite = false } = req.body;
      
      if (!Array.isArray(smiles) || smiles.length === 0) {
        return res.status(400).json({ 
          error: 'SMILES array required' 
        });
      }
      
      logger.info(`SDF generation request for ${smiles.length} molecules`);
      
      const result = await molecularProcessor.processSmiles(smiles, overwrite);
      res.json(result);
    } catch (error) {
      next(error);
    }
  });
  
  app.post('/api/name-to-sdf', async (req, res, next) => {
    try {
      const { name, overwrite = false } = req.body;

      if (!name || typeof name !== 'string') {
        return res.status(400).json({
          error: 'Chemical name required'
        });
      }

      logger.info(`Name to SDF request: "${name}"`);

      const result = await molecularProcessor.generateSDFByName(name, overwrite);

      if (!result) {
        return res.status(404).json({
          error: `Could not find chemical: ${name}`
        });
      }

      res.json(result);
    } catch (error) {
      next(error);
    }
  });

  // Test data endpoint
  app.get('/api/test-data/:testName', async (req, res, next) => {
    try {
      const { testName } = req.params;
      const testDataPath = path.join(__dirname, '..', '..', '..', 'tests', 'fixtures', 'visual-test-data', `${testName}.json`);

      if (!fs.existsSync(testDataPath)) {
        return res.status(404).json({
          error: `Test data not found: ${testName}`
        });
      }

      const testData = JSON.parse(fs.readFileSync(testDataPath, 'utf8'));
      res.json(testData);
    } catch (error) {
      next(error);
    }
  });
}


async function startServer(container) {
  const config = await container.get('config');
  const logger = await container.get('logger');
  

  try {
    // Validate configuration (silent on success)
    config.validate();

    // CRITICAL: Validate OPENAI_API_KEY for local development
    if (!config.isProduction() && !config.get('openai.apiKey')) {
      logger.error('❌ OPENAI_API_KEY not found in environment');
      logger.error('   Create .env file in project root with:');
      logger.error('   OPENAI_API_KEY=sk-your-key-here');
      throw new Error('OPENAI_API_KEY required for local development - check .env file');
    }

    // Initialize database if enabled
    if (config.get('database.enabled')) {
      const database = await container.get('database');
      if (database && database.initialize) {
        await database.initialize();
      }
    }
    
    // Create Express app
    const app = await createApp(container);
    
    // Start HTTP server
    const port = config.get('port') || 8080;
    const server = app.listen(port, () => {
      logger.info('Server configuration:', {
        environment: config.get('nodeEnv'),
        databaseEnabled: config.get('database.enabled'),
        port
      });
    });
    
    // Start HTTPS server if configured
    if (config.get('ssl.certPath') && config.get('ssl.keyPath')) {
      const HttpsServer = require('./https-server');
      const httpsServer = new HttpsServer(app);
      const httpsPort = config.get('ssl.httpsPort') || 3001;
      
      httpsServer.start(httpsPort);
      logger.info(`🔒 HTTPS server running on port ${httpsPort}`);
    }
    
    // Graceful shutdown
    process.on('SIGTERM', async () => {
      logger.info('SIGTERM received, shutting down gracefully...');
      
      server.close(() => {
        logger.info('HTTP server closed');
      });
      
      // Cleanup services
      const database = await container.get('database');
      if (database && database.close) {
        await database.close();
        logger.info('Database connections closed');
      }
      
      process.exit(0);
    });
    
    return server;
  } catch (error) {
    logger.error('❌ Server startup failed:', error);
    process.exit(1);
  }
}

// If this file is run directly, start the server
if (require.main === module) {
  const container = createContainer();
  startServer(container).catch(console.error);
}

// Export for testing
module.exports = {
  createApp,
  startServer,
  setupMolecularRoutes: setupChemicalPredictionRoutes,
};
