const net = require("net");
const path = require("path");

const findAvailablePort = async (startPort) => {
  const isPortAvailable = (port) => {
    return new Promise((resolve) => {
      const tester = net
        .createServer()
        .once("error", () => resolve(false))
        .once("listening", () => {
          tester.once("close", () => resolve(true)).close();
        })
        .listen(port);
    });
  };

  let port = startPort;
  while (port < startPort + 100) {
    if (await isPortAvailable(port)) {
      return port;
    }
    port++;
  }
  throw new Error(
    `No available ports found between ${startPort} and ${startPort + 100}`,
  );
};

const attemptPortCleanup = async (port) => {
  try {
    const { execSync } = require("child_process");

    // Try to find and kill processes using the port (Unix/macOS)
    if (process.platform !== "win32") {
      try {
        const result = execSync(`lsof -ti:${port}`, { encoding: "utf8", stdio: "pipe" });
        const pids = result.trim().split("\n").filter(Boolean);

        if (pids.length > 0) {
          console.log(`🧹 Found ${pids.length} process(es) using port ${port}, attempting cleanup...`);

          for (const pid of pids) {
            try {
              execSync(`kill -9 ${pid}`, { stdio: "pipe" });
              console.log(`   ✅ Killed process ${pid}`);
            } catch (e) {
              console.log(`   ⚠️ Could not kill process ${pid} (may not have permission)`);
            }
          }

          // Wait a moment for cleanup
          await new Promise(resolve => setTimeout(resolve, 500));
          return true;
        }
      } catch (e) {
        // No processes found or lsof not available
        return false;
      }
    }
    return false;
  } catch (error) {
    console.log(`⚠️ Port cleanup failed: ${error.message}`);
    return false;
  }
};

const cleanupPorts = async () => {
  const portsToCheck = [8080, 3003, 3004, 3005];

  for (const port of portsToCheck) {
    try {
      const server = net.createServer();
      const isPortFree = await new Promise((resolve) => {
        server.listen(port, "127.0.0.1", () => {
          server.close();
          resolve(true);
        });
        server.on("error", () => resolve(false));
      });

      if (!isPortFree) {
        console.log(`⚠️ Port ${port} is in use, attempting cleanup...`);
        // Kill any processes using this port
        try {
          const { execSync } = require("child_process");
          execSync(`lsof -ti:${port} | xargs kill -9 2>/dev/null || true`);
          console.log(`✅ Cleaned up port ${port}`);
        } catch (e) {
          console.log(`⚠️ Could not cleanup port ${port}: ${e.message}`);
        }
      }
    } catch (e) {
      // Ignore cleanup errors
    }
  }
};

const getLocalIPAddress = () => {
  const os = require("os");
  const interfaces = os.networkInterfaces();
  for (const name of Object.keys(interfaces)) {
    for (const iface of interfaces[name]) {
      if (iface.family === "IPv4" && !iface.internal) {
        return iface.address;
      }
    }
  }
  return "127.0.0.1";
};

const startServers = async (app, configuration, initializeDatabase) => {
  const isTestMode =
    configuration.get('nodeEnv') === "test" || !!process.env.JEST_WORKER_ID;
  const isIntegrationTest = !!process.env.INTEGRATION_TEST;

  // Store server instances for cleanup
  let httpServer;
  let httpsServerInstance;

  if (!isTestMode || isIntegrationTest) {

    const startServer = async () => {
      await cleanupPorts();

      try {
        let actualPort = configuration.get('port');

        if (actualPort === 8080) {
          try {
            if (!(await (async () => {
              const net = require("net");
              return new Promise((resolve) => {
                const server = net.createServer();
                server.listen(actualPort, "0.0.0.0", () => {
                  server.once("close", () => resolve(true));
                  server.close();
                });
                server.on("error", () => resolve(false));
              });
            })())) {
              console.log(`⚠️ Port ${actualPort} is in use, attempting cleanup...`);
              const cleanupSuccessful = await attemptPortCleanup(actualPort);

              if (cleanupSuccessful) {
                console.log(`✅ Port ${actualPort} cleanup completed, retrying...`);
                await new Promise(resolve => setTimeout(resolve, 1000));
              }
            }

            actualPort = await findAvailablePort(actualPort);
            if (actualPort !== configuration.get('port')) {
              console.log(
                `⚠️ Port ${configuration.get('port')} is in use, using port ${actualPort} instead`,
              );
            }
          } catch (error) {
            console.error(`❌ Could not find available port: ${error.message}`);
            console.log(`💡 Try: pkill -f "node.*server.js" or use a different port`);
            process.exit(1);
          }
        }

        const localIP = getLocalIPAddress();
        const endpoints = [
          `http://localhost:${actualPort}`,
          ...(localIP && localIP !== '127.0.0.1' ? [`http://${localIP}:${actualPort}`] : [])
        ];

        httpServer = app.listen(actualPort, "0.0.0.0", () => {
          const logger = require('./file-logger');
          endpoints.forEach((endpoint) => {
            logger.startup(`HTTP available at ${endpoint}`);
          });
          
          // Auto-inject test case if specified
          const autoTestInject = process.env.AUTO_TEST_INJECT;
          if (autoTestInject) {
            logger.info(`🧪 Auto-injection enabled for test case: ${autoTestInject}`);
            // Wait a moment for server to be fully ready, then inject
            setTimeout(async () => {
              try {
                const { injectPrediction } = require('../../../scripts/inject-test');
                const testCases = require('../../../scripts/inject-test').TEST_CASES;
                
                if (testCases[autoTestInject]) {
                  logger.info(`🔥 Auto-injecting test: ${autoTestInject}`);
                  await injectPrediction(testCases[autoTestInject]);
                  logger.success(`✅ Auto-injection completed for: ${autoTestInject}`);
                } else {
                  logger.warn(`⚠️ Unknown test case for auto-injection: ${autoTestInject}`);
                  logger.info(`Available: ${Object.keys(testCases).join(', ')}`);
                }
              } catch (error) {
                logger.warn(`⚠️ Auto-injection failed: ${error.message}`);
              }
            }, 2000); // 2 second delay to ensure server is ready
          }
        });
      } catch (error) {
        console.error(`❌ Failed to start server: ${error.message}`);
        process.exit(1);
      }
    };

    const startHttpsServer = async () => {
      try {
        const HttpsServer = require("../api/https-server");
        const httpsServer = new HttpsServer(app, configuration.get('ssl.httpsPort'));
        const result = await httpsServer.start();
        httpsServerInstance = result ? result.server : null;
        const httpsEndpoints = result ? result.endpoints : [];

        if (httpsServerInstance) {
          httpsServerInstance.on("error", (error) => {
            const logger = require('./file-logger');
            logger.error("HTTPS server error after startup", { error: error.message });
          });

          try {
            const cleanupRegistry = require('../../tests/support/fixtures/cleanup-registry');
            cleanupRegistry.register(httpsServerInstance);
          } catch (e) {
            // Cleanup registry not available
          }

          httpsEndpoints.forEach((endpoint) => {
            const logger = require('./file-logger');
            logger.startup(`HTTPS available at ${endpoint}`);
          });
        } else {
          console.log("⚠️ HTTPS server not started - continuing with HTTP only");
        }
      } catch (error) {
        console.error("❌ Failed to start HTTPS server:", error.message);
        console.log("💡 Continuing with HTTP server only");
      }
    };

    await Promise.all([
      startHttpsServer(),
      startServer()
    ]);

    if (httpServer) {
      httpServer.on("error", (error) => {
        if (error.code === "EADDRINUSE") {
          console.error(`❌ Port ${configuration.get('port')} is already in use`);
          console.log(`💡 Solutions:`);
          console.log(`   1. Kill existing process: pkill -f "node.*server.js"`);
          console.log(`   2. Use different port: PORT=8081 npm start`);
          console.log(`   3. Check what's using the port: lsof -i :${configuration.get('port')}`);
          process.exit(1);
        } else if (error.code === "EACCES") {
          console.error(`❌ Permission denied: Cannot bind to port ${configuration.get('port')}`);
          console.log(`💡 Try using a port > 1024 or run with sudo`);
          process.exit(1);
        } else {
          console.error("❌ Server error:", error.message);
          console.log(`💡 Check your network configuration and try again`);
          process.exit(1);
        }
      });
    }

    // Initialize database after servers start
    if (initializeDatabase) {
      await initializeDatabase();
    }
  }

  // Graceful shutdown handling
  if (!isTestMode) {
    const gracefulShutdown = (signal) => {
      const closeServer = (server, name) => {
        return new Promise((resolve) => {
          if (server) {
            server.close(() => {
              resolve();
            });
          } else {
            resolve();
          }
        });
      };

      Promise.all([
        closeServer(httpServer, "HTTP"),
        closeServer(httpsServerInstance, "HTTPS"),
      ]).then(() => {
        process.exit(0);
      });

      setTimeout(() => {
        process.exit(1);
      }, 5000);
    };

    process.on("SIGTERM", () => gracefulShutdown("SIGTERM"));
    process.on("SIGINT", () => gracefulShutdown("SIGINT"));
    process.on("SIGUSR2", () => gracefulShutdown("SIGUSR2"));
  }

  return { httpServer, httpsServer: httpsServerInstance };
};

module.exports = {
  findAvailablePort,
  attemptPortCleanup,
  cleanupPorts,
  getLocalIPAddress,
  startServers,
};