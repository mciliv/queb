const { createContainer } = require('../../core/services');
const { createApp, setupChemicalPredictionRoutes } = require('./app');

// Main server startup - consolidates app creation + server lifecycle
async function startServer(container) {
  // Resolve core services once
  const config = await container.get('config');
  const logger = await container.get('logger');

  try {
    // Create the Express app
    const app = await createApp({ config, logger, container });

    // Start HTTP server
    const port = config.get('port') || 8080;
    const server = app.listen(port, () => {
      logger.info('Server configuration:', {
        environment: config.get('nodeEnv'),
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

// Export for testing - createApp remains available for isolated testing
module.exports = {
  createApp,           // For testing app creation without server
  startServer,         // For testing full server lifecycle
  setupMolecularRoutes: setupChemicalPredictionRoutes,
};
