const https = require("https");
const fs = require("fs");
const path = require("path");
const os = require("os");
const { execSync } = require("child_process");
const net = require("net");
const WebSocket = require("ws");

// Simple logging utility for consistency
const log = {
  success: (msg) => console.log(msg),
  warning: (msg) => console.warn(msg),
  error: (msg) => console.error(msg)
};

// Register with cleanup system if available
let cleanupRegistry = null;
try {
  cleanupRegistry = require('../../test/fixtures/cleanup-registry');
} catch (e) {
  // Cleanup registry not available in production
}

class HttpsServer {
  constructor(app, port = 3001) {
    this.app = app;
    this.requestedPort = port;
    this.actualPort = port;
    this.localIP = this.getLocalIPAddress();
  }

  getLocalIPAddress() {
    const interfaces = os.networkInterfaces();
    for (const name of Object.keys(interfaces)) {
      for (const iface of interfaces[name]) {
        if (iface.family === "IPv4" && !iface.internal) {
          return iface.address;
        }
      }
    }
    return "127.0.0.1";
  }

  // Check if a port is available
  async isPortAvailable(port) {
    return new Promise((resolve) => {
      const server = net.createServer();
      
      server.listen(port, "0.0.0.0", () => {
        server.once("close", () => resolve(true));
        server.close();
      });
      
      server.on("error", () => resolve(false));
    });
  }

  // Find an available port - NO BACKUP PORTS (only try the requested port)
  async findAvailablePort(startPort) {
    if (await this.isPortAvailable(startPort)) {
      return startPort;
    }
    throw new Error(`Port ${startPort} is not available - no backup ports allowed`);
  }

  loadCertificates() {
    // Try multiple certificate locations for flexibility
    const certLocations = [
      path.join(__dirname, "certs"),                    // backend/api/certs (legacy)
      path.join(process.cwd(), "certs"),               // project root certs
      path.join(process.cwd(), "config", "certs")      // config/certs
    ];

    for (const certDir of certLocations) {
      const keyPath = path.join(certDir, "key.pem");
      const certPath = path.join(certDir, "cert.pem");

      if (fs.existsSync(keyPath) && fs.existsSync(certPath)) {
        try {
          log.success(`✅ Using mkcert SSL certificates from: ${certDir}`);
          return { key: fs.readFileSync(keyPath), cert: fs.readFileSync(certPath) };
        } catch (error) {
          log.warning(`⚠️ Could not read certificates from ${certDir}: ${error.message}`);
          continue;
        }
      }
    }

    // No certificates found
    log.error("❌ No SSL certificates found");
    log.warning("💡 To generate trusted certificates, run:");
    log.warning("   node ../../web/generate-certs.js");

    return null;
  }

  async start() {
    try {
      const credentials = this.loadCertificates();
      if (!credentials) {
        log.warning("⚠️ HTTPS not started — SSL certificates not found");
        return null;
      }

      // Find an available port
      try {
        this.actualPort = await this.findAvailablePort(this.requestedPort);
        
        if (this.actualPort !== this.requestedPort) {
          log.warning(
            `⚠️ HTTPS port ${this.requestedPort} in use, using port ${this.actualPort} instead`
          );
        }
      } catch (error) {
        log.error(`❌ Could not find available HTTPS port: ${error.message}`);
        log.warning("💡 Try stopping other services or use a different port range");
        return null;
      }

      const server = https.createServer(credentials, this.app);

      // Add WebSocket support for browser extensions/tools
      const wss = new WebSocket.Server({ server });

      wss.on('connection', (ws) => {
        // Handle WebSocket connections (basic implementation)
        ws.on('message', (message) => {
          // Echo back any messages received
          ws.send('WebSocket connection established');
        });

        ws.on('error', (error) => {
          // Silently handle WebSocket errors
          console.log('WebSocket connection error:', error.message);
        });
      });

      // Register with cleanup system if available
      if (cleanupRegistry) {
        cleanupRegistry.register(server);
      }
      
      return new Promise((resolve, reject) => {
        server.listen(this.actualPort, "0.0.0.0", () => {
          console.log(`https://localhost:${this.actualPort}`);
          if (this.localIP && this.localIP !== "127.0.0.1") {
            console.log(`https://${this.localIP}:${this.actualPort}`);
          }
          resolve(server);
        });

        server.on("error", (error) => {
          if (error.code === "EADDRINUSE") {
            log.error(`❌ HTTPS port ${this.actualPort} is already in use`);
            log.warning("💡 No backup ports allowed - stopping server");
            server.close();
            resolve(null);
          } else if (error.code === "EACCES") {
            log.error(`❌ Permission denied: Cannot bind to HTTPS port ${this.actualPort}`);
            log.warning("💡 Try using a port > 1024 or run with appropriate permissions");
            resolve(null);
          } else {
            log.error("❌ HTTPS server error:", error.message);
            log.warning("💡 HTTPS server will continue without SSL");
            resolve(null);
          }
        });
      });

    } catch (error) {
      log.error("❌ Failed to start HTTPS server:", error.message);
      log.warning("💡 Continuing without HTTPS support");
      return null;
    }
  }

  // Graceful shutdown
  async stop(server) {
    if (server) {
      return new Promise((resolve) => {
        server.close(() => {
          log.success("🔒 HTTPS server stopped gracefully");
          resolve();
        });
        
        // Force close after timeout
        setTimeout(() => {
          if (server.listening) {
            server.destroy();
            log.warning("⚠️ HTTPS server force closed after timeout");
          }
          resolve();
        }, 5000);
      });
    }
  }
}

module.exports = HttpsServer;
