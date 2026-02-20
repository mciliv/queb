const express = require('express');
const cors = require('cors');
const path = require('path');
const fs = require('fs');
const createError = require('http-errors');

/**
 * Create and configure Express application with injected services.
 * This file contains the "general app logic": middleware, routes, SPA/static serving, and error handling.
 *
 * @param {Object} services - Resolved core services
 * @param {Object} services.config - Configuration service
 * @param {Object} services.logger - Logger service
 * @param {ServiceContainer} services.container - DI container for other services
 * @returns {Promise<express.Application>} Configured Express app
 */
async function createApp({ config, logger, container }) {
  const app = express();

  // Get remaining services
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
        .replace(/\{\{OG_IMAGE\}\}/g, '/assets/favicon.svg')
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

    res.send = function (data) {
      res.send = originalSend;
      const duration = Date.now() - start;
      logger.info(`${req.method} ${req.path} - ${res.statusCode} (${duration}ms)`);
      return res.send(data);
    };

    next();
  });

  // Health check (support both /api/health and /health for compatibility)
  const healthHandler = (req, res) => {
    res.json({
      status: 'ok',
      environment: config.get('nodeEnv'),
      version: process.env.npm_package_version || '1.0.0',
    });
  };
  app.get('/api/health', healthHandler);

  app.get('/health', healthHandler);

  // Frontend error logging endpoint (dev only)
  // THIN CONTROLLER: HTTP ↔ Domain translation only. No domain logic, no infra coupling.
  // NOTE: Use Cursor's builtin automatic screenshots for error debugging and reproduction
  app.post('/api/log-error', async (req, res) => {
    try {
      // ONE DOMAIN CALL: All domain logic abstracted into use-case
      await container.get('logErrorUseCase').execute(req.body);
      res.status(200).json({ received: true });
    } catch (error) {
      // HTTP error handling only - domain errors translated here
      res.status(500).json({ error: 'Failed to log error' });
    }
  });

  await setupChemicalPredictionRoutes(app, { logger, container });
  await setupDatabaseRecommendationRoutes(app, { logger, container });

  // Setup hotel routes (secure video hosting and request card)
  const { setupHotelRoutes } = require('../routes/hotel');
  setupHotelRoutes(app, { logger });

  // Setup chemical analysis routes (AI-powered compound analysis)
  const { setupChemicalAnalysisRoutes } = require('../routes/chemical-analysis');
  setupChemicalAnalysisRoutes(app, { logger, container });

  // Setup user routes if database is enabled
  if (config.get('database.enabled')) {
    const { setupUserRoutes } = require('./user-routes');
    await setupUserRoutes(app, { config, logger, container });
  }

  // SPA catch-all - serve index.html for non-API, non-file requests
  app.get('*', (req, res, next) => {
    // Skip if it's an API route, hotel route, or file request
    if (req.path.startsWith('/api/') || req.path.startsWith('/hotel/') ||
        req.path.startsWith('/sdf_files/') || req.path.includes('.')) {
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
    const status = err.status || err.statusCode || 500;

    // Don't log 404s for common browser asset requests
    if (status === 404) {
      const isAssetRequest = /^\/(favicon\.ico|apple-touch-icon.*\.png|images\/favicon\.svg)$/.test(req.path);
      if (isAssetRequest) {
        return res.sendStatus(404);
      }
    }
    
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
        context: err.context,
      };
    } else {
      handled = await errorHandler.handle(err, {
        path: req.path,
        method: req.method,
        ip: req.ip,
      });
    }


    // For API routes, always return JSON
    if (req.path.startsWith('/api/')) {
      return res.status(status).json({
        error: handled.message || err.message || 'Internal server error',
        code: handled.code,
        recoverable: handled.recoverable,
        recovery: handled.recovery,
        timestamp: handled.timestamp,
        ...(config.isDevelopment() && {
          stack: err.stack,
          context: handled.context,
          originalMessage: err.message,
        }),
      });
    }

    // For non-API routes, return appropriate response based on status
    if (status === 404) {
      return res.status(404).send('Not found');
    }

    res.status(status).json({
      error: handled.message || err.message || 'Internal server error',
      ...(config.isDevelopment() && {
        stack: err.stack,
        context: handled.context,
      }),
    });
  });

  return app;
}

async function setupChemicalPredictionRoutes(app, { logger, container }) {
  const structuralizer = await container.get('structuralizer');
  const molecularProcessor = await container.get('molecularProcessor');

  app.post('/api/structuralize', async (req, res, next) => {
    try {
      const { text, lookupMode = 'GPT-5' } = req.body;

      if (!text || typeof text !== 'string') {
        return res.status(400).json({
          error: 'Missing or invalid "text" parameter',
        });
      }

      if (typeof lookupMode !== 'string') {
        return res.status(400).json({
          error: 'Invalid "lookupMode" parameter',
        });
      }

      logger.info(`Text prediction request: "${text}" (mode: ${lookupMode})`);

      const result = await structuralizer.chemicals({
        object: text.trim(),
        lookupMode,
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
          error: 'Missing imageBase64 parameter',
        });
      }

      logger.info('Image prediction request', {
        hasCoordinates: x !== undefined && y !== undefined,
      });

      const result = await structuralizer.chemicals({
        imageBase64,
        x,
        y,
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
          error: 'SMILES array required',
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
          error: 'Chemical name required',
        });
      }

      logger.info(`Name to SDF request: "${name}"`);

      const result = await molecularProcessor.generateSDFByName(name, overwrite);

      if (!result) {
        return res.status(404).json({
          error: `Could not find chemical: ${name}`,
        });
      }

      res.json(result);
    } catch (error) {
      next(error);
    }
  });

  // 3D structure endpoint - returns SDF data for AI display
  app.get('/api/3d-structure/:name', async (req, res, next) => {
    try {
      const { name } = req.params;
      const recordType = req.query.record_type || '3d';

      if (!name || typeof name !== 'string') {
        return res.status(400).json({
          error: 'Chemical name is required',
        });
      }

      logger.info(`3D structure request: "${name}" (type: ${recordType})`);

      // Fetch SDF from PubChem
      const encoded = encodeURIComponent(name);
      const pubchemUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encoded}/SDF?record_type=${recordType}`;
      
      const fetch = (await import('node-fetch')).default;
      const response = await fetch(pubchemUrl);

      if (!response.ok) {
        if (response.status === 404) {
          return res.status(404).json({
            error: `3D structure not found for: ${name}`,
            suggestion: 'Try common chemical names like "caffeine", "aspirin", or "glucose"',
          });
        }
        throw new Error(`PubChem API error: ${response.status} ${response.statusText}`);
      }

      const sdfContent = await response.text();

      // Get compound metadata for context
      const metadataUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encoded}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IUPACName/JSON`;
      let metadata = null;
      try {
        const metaResponse = await fetch(metadataUrl);
        if (metaResponse.ok) {
          const metaData = await metaResponse.json();
          const props = metaData?.PropertyTable?.Properties?.[0];
          if (props) {
            metadata = {
              cid: props.CID,
              formula: props.MolecularFormula,
              molecularWeight: props.MolecularWeight,
              smiles: props.CanonicalSMILES,
              iupacName: props.IUPACName,
            };
          }
        }
      } catch (_) {
        // Metadata fetch is optional, continue without it
      }

      // Return SDF content with metadata
      res.set('Content-Type', 'application/json');
      res.json({
        name: name,
        format: 'sdf',
        record_type: recordType,
        sdf_content: sdfContent,
        metadata: metadata,
        note: 'Use the sdf_content field to display the 3D structure. The SDF format is compatible with most molecular viewers.',
      });
    } catch (error) {
      next(error);
    }
  });

  // Test data endpoint
  app.get('/api/test-data/:testName', async (req, res, next) => {
    try {
      const { testName } = req.params;
      const testDataPath = path.join(
        __dirname,
        '..',
        '..',
        '..',
        'tests',
        'fixtures',
        'visual-test-data',
        `${testName}.json`,
      );

      if (!fs.existsSync(testDataPath)) {
        return res.status(404).json({
          error: `Test data not found: ${testName}`,
        });
      }

      const testData = JSON.parse(fs.readFileSync(testDataPath, 'utf8'));
      res.json(testData);
    } catch (error) {
      next(error);
    }
  });
}

async function setupDatabaseRecommendationRoutes(app, { logger, container }) {
  const DatabaseRecommender = require('../services/database-recommender');
  const aiService = await container.get('aiService');
  const nameResolver = require('../services/name-resolver');

  const recommender = new DatabaseRecommender({
    logger,
    aiService,
    nameResolver,
  });

  /**
   * Use AI to select best database and query for chemical contents
   * POST /api/find-chemical-contents
   * Body: { item: "coffee" }
   * Returns: { item, selectedDatabase, chemicals, metadata, queryInfo }
   */
  app.post('/api/find-chemical-contents', async (req, res, next) => {
    try {
      const { item } = req.body;

      if (!item || typeof item !== 'string') {
        return res.status(400).json({
          error: 'Missing or invalid "item" parameter',
        });
      }

      logger.info(`Chemical contents search request: "${item}"`);

      const result = await recommender.findChemicalContents(item);

      res.json(result);
    } catch (error) {
      next(error);
    }
  });

  /**
   * Recommend optimal database(s) for searching contents of an item (legacy endpoint)
   * POST /api/recommend-database
   * Body: { item: "coffee" }
   * Returns: { item, analysis, recommendations, databases }
   */
  app.post('/api/recommend-database', async (req, res, next) => {
    try {
      const { item } = req.body;

      if (!item || typeof item !== 'string') {
        return res.status(400).json({
          error: 'Missing or invalid "item" parameter',
        });
      }

      logger.info(`Database recommendation request: "${item}"`);

      const result = await recommender.recommendDatabase(item);

      res.json(result);
    } catch (error) {
      next(error);
    }
  });
}

module.exports = {
  createApp,
  setupChemicalPredictionRoutes,
  setupDatabaseRecommendationRoutes,
};



