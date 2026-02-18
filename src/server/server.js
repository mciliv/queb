const { createContainer } = require('../../core/ServiceContainer');
const { createApp, setupChemicalPredictionRoutes } = require('./app');

async function startServer(container) {
  const config = await container.get('config');
  const logger = await container.get('logger');

  try {
    const app = await createApp({ config, logger, container });
    const server = httpServers(app);

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

  function httpServers(app) {
    const port = config.get('port') || 8080;
    const server = app.listen(port, () => {
      logger.info('Server configuration:', {
        environment: config.get('nodeEnv'),
        port
      });
    });

    if (config.get('ssl.certPath') && config.get('ssl.keyPath')) {
      const HttpsServer = require('./https-server');
      const httpsServer = new HttpsServer(app);
      const httpsPort = config.get('ssl.httpsPort') || 3001;

      httpsServer.start(httpsPort);
      logger.info(`🔒 HTTPS server running on port ${httpsPort}`);
    }
    return server;
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
