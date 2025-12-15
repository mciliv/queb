const { createContainer } = require('../../core/services');
const { createApp, setupChemicalPredictionRoutes } = require('./app');
const { validateLocalDevEnv } = require('./validate-env');

async function startServer(container) {
  const config = await container.get('config');
  const logger = await container.get('logger');
  

  try {
    // Validate configuration (silent on success)
    config.validate();

    // Environment validations that are specific to local/dev startup
    validateLocalDevEnv(config, logger);

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
