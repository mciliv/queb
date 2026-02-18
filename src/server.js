#!/usr/bin/env node

const express = require('express');
const cors = require('cors');
const path = require('path');
const fs = require('fs');

// Load environment variables
require('dotenv').config();

// Import essential services
const AIService = require('./src/server/services/AIService');
const Structuralizer = require('./src/server/services/Structuralizer');
const MolecularProcessor = require('./src/server/services/molecular-processor');
const { resolveName } = require('./src/server/services/name-resolver');
const promptEngine = require('./src/core/PromptEngine'); // Still needed by Structuralizer
const SimpleCache = require('./src/server/services/SimpleCache');
const { fetchHighestQualityIngredients } = require('./src/server/services/ingredient-delivery');

// CLI argument parsing
const args = process.argv.slice(2);
const parsedArgs = {};
for (let i = 0; i < args.length; i++) {
  const arg = args[i];
  if (arg.startsWith('--')) {
    const [key, value] = arg.slice(2).split('=');
    parsedArgs[key] = value || true;
  }
}

const PORT = parsedArgs.port ? parseInt(parsedArgs.port) : (process.env.PORT || 8080);
const testMode = parsedArgs.test;
const devMode = parsedArgs.dev;
global.serverConfig = { testMode, devMode, port: PORT };

// Initialize services
const aiService = new AIService();
const molecularProcessor = new MolecularProcessor();

const logger = {
  info: (...args) => console.log('[INFO]', ...args),
  error: (...args) => console.error('[ERROR]', ...args),
  warn: (...args) => console.warn('[WARN]', ...args)
};

const errorHandler = {
  handle: async (err, context) => {
    logger.error('Error:', err.message, context);
    return {
      message: err.message || 'Internal server error',
      code: err.code || 'INTERNAL_ERROR',
      recoverable: false,
      recovery: null,
      timestamp: new Date().toISOString(),
      context
    };
  }
};

// User prompt: "extract the cache components into it's own file"
// Initialize simple in-memory cache
const cache = new SimpleCache({
  maxSize: 1000,          // Max 1000 cached entries
  defaultTTL: 300000,     // 5 minutes default
  cleanupInterval: 60000  // Cleanup expired entries every minute
});

// User prompt: "inject it?" - Yes! Create and inject CacheWrapper
const CacheWrapper = require('./src/server/services/CacheWrapper');
const cacheWrapper = new CacheWrapper({
  cache: cache,
  enabled: true,
  logger: logger,
  defaultTTL: 300000
});

const structuralizer = new Structuralizer({
  aiService,
  molecularProcessor,
  nameResolver: { resolveName },
  promptEngine,
  errorHandler,
  logger,
  cache: cache,  // Still available for backward compatibility
  cacheWrapper: cacheWrapper,  // Injected CacheWrapper instance
  config: {
    aiTimeout: 10000,
    maxRetries: 2,
    cacheEnabled: true  // Enable caching
  }
});

const app = express();

// Middleware
app.use(cors());
app.use(express.json({ limit: '50mb' }));
app.use(express.urlencoded({ extended: true, limit: '50mb' }));

// Serve static files
const distPath = path.join(__dirname, 'src', 'client', 'dist');
const publicPath = path.join(__dirname, 'public');
const assetsPath = path.join(__dirname, 'src', 'client', 'assets');

app.use(express.static(publicPath));
app.use('/assets', express.static(assetsPath));
app.use('/dist', express.static(distPath));
app.use(express.static(distPath, { index: false }));

// Serve SDF files
const sdfDir = path.join(__dirname, 'src', 'public', 'sdf_files');
if (fs.existsSync(sdfDir)) {
  app.use('/sdf_files', express.static(sdfDir));
}

// Health check
app.get('/api/health', (req, res) => {
  res.json({
    status: 'ok',
    timestamp: new Date().toISOString(),
    config: global.serverConfig,
    cache: cache.getStats()
  });
});

// Cache statistics endpoint
app.get('/api/cache/stats', (req, res) => {
  res.json({
    stats: cache.getStats(),
    enabled: true,
    type: 'simple-in-memory'
  });
});

// Clear cache endpoint (for development)
app.post('/api/cache/clear', async (req, res) => {
  await cache.clear();
  logger.info('Cache cleared manually');
  res.json({ 
    success: true, 
    message: 'Cache cleared',
    stats: cache.getStats()
  });
});

app.post('/api/ingredients/highest-quality', async (req, res) => {
  try {
    const { city, lat, lng, limit } = req.body || {};

    if (city && typeof city !== 'string') {
      return res.status(400).json({ error: 'Invalid "city" parameter' });
    }

    const hasLatLng = lat !== undefined || lng !== undefined;
    if (hasLatLng && (typeof lat !== 'number' || typeof lng !== 'number')) {
      return res.status(400).json({ error: 'Invalid "lat"/"lng" parameters' });
    }

    if (!city && !hasLatLng) {
      return res.status(400).json({ error: 'Provide "city" or "lat"/"lng"' });
    }

    if (limit !== undefined && (!Number.isInteger(limit) || limit <= 0)) {
      return res.status(400).json({ error: 'Invalid "limit" parameter' });
    }

    const results = await fetchHighestQualityIngredients({
      city: city ? city.trim() : undefined,
      lat,
      lng,
      limit
    });

    res.json({ results });
  } catch (error) {
    logger.error('Ingredient delivery error:', error);
    res.status(500).json({ error: error.message || 'Ingredient delivery failed' });
  }
});

app.post('/api/structuralize', async (req, res) => {
  try {
    const { text, lookupMode = 'GPT-5' } = req.body;

    if (!text || typeof text !== 'string') {
      return res.status(400).json({
        error: 'Missing or invalid "text" parameter'
      });
    }

    logger.info(`Text analysis request: "${text}" (mode: ${lookupMode})`);

    const result = await structuralizer.chemicals({
      object: text.trim(),
      lookupMode
    });

    res.json(result);
  } catch (error) {
    logger.error('Structuralize error:', error);
    res.status(500).json({
      error: error.message || 'Analysis failed',
      details: process.env.NODE_ENV === 'development' ? error.stack : undefined
    });
  }
});

app.post('/api/structuralize-image', async (req, res) => {
  try {
    const { imageBase64, x, y } = req.body;

    if (!imageBase64) {
      return res.status(400).json({
        error: 'Missing imageBase64 parameter'
      });
    }

    logger.info('Image analysis request', {
      hasCoordinates: x !== undefined && y !== undefined
    });

    const result = await structuralizer.chemicals({
      imageBase64,
      x,
      y
    });

    res.json(result);
  } catch (error) {
    logger.error('Image analysis error:', error);
    res.status(500).json({
      error: error.message || 'Image analysis failed',
      details: process.env.NODE_ENV === 'development' ? error.stack : undefined
    });
  }
});

app.post('/api/generate-sdfs', async (req, res) => {
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
    logger.error('SDF generation error:', error);
    res.status(500).json({
      error: error.message || 'SDF generation failed',
      details: process.env.NODE_ENV === 'development' ? error.stack : undefined
    });
  }
});

// POST /api/analyze - Simple text analysis (legacy/alternative endpoint)
const { setupChemicalAnalysisRoutes } = require('./src/server/routes/chemical-analysis');
setupChemicalAnalysisRoutes(app, { logger, container: null });

// Serve index.html
app.get('/', (req, res) => {
  try {
    const indexPath = path.join(distPath, 'index.html');
    if (fs.existsSync(indexPath)) {
      let htmlContent = fs.readFileSync(indexPath, 'utf8');
      const cacheBust = Date.now();
      htmlContent = htmlContent.split('{{CACHE_BUST}}').join(cacheBust);
      res.type('html').send(htmlContent);
    } else {
      res.status(500).send('Frontend not built. Run: npm run build');
    }
  } catch (error) {
    logger.error('Failed to serve index.html:', error);
    res.status(500).send('Server error');
  }
});

// SPA catch-all
app.get('*', (req, res) => {
  if (req.path.startsWith('/api/') || req.path.startsWith('/sdf_files/') || req.path.includes('.')) {
    return res.status(404).json({ error: 'Not found' });
  }

  try {
    const indexPath = path.join(distPath, 'index.html');
    if (fs.existsSync(indexPath)) {
      let htmlContent = fs.readFileSync(indexPath, 'utf8');
      const cacheBust = Date.now();
      htmlContent = htmlContent.split('{{CACHE_BUST}}').join(cacheBust);
      res.type('html').send(htmlContent);
    } else {
      res.status(500).send('Frontend not built. Run: npm run build');
    }
  } catch (error) {
    logger.error('Failed to serve index.html:', error);
    res.status(500).send('Server error');
  }
});

// Start server
app.listen(PORT, () => {
  console.log(`🚀 Chemical Analyzer running on http://localhost:${PORT}`);
  console.log(`📝 Enter text like "coffee", "water", "aspirin" to analyze chemical content`);
});
