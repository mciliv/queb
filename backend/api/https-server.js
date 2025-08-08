const https = require("https");
const fs = require("fs");
const path = require("path");
const os = require("os");
const { execSync } = require("child_process");
const net = require("net");

// Simple logging utility for consistency
const log = {
  success: (msg) => console.log(msg),
  warning: (msg) => console.log(msg),
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

  // Find an available port starting from the requested port
  async findAvailablePort(startPort) {
    for (let port = startPort; port < startPort + 100; port++) {
      if (await this.isPortAvailable(port)) {
        return port;
      }
    }
    throw new Error(`No available ports found in range ${startPort}-${startPort + 99}`);
  }

  loadCertificates() {
    const certDir = path.join(__dirname, "certs");
    const keyPath = path.join(certDir, "key.pem");
    const certPath = path.join(certDir, "cert.pem");

    // Check for mkcert certificates
    if (fs.existsSync(keyPath) && fs.existsSync(certPath)) {
      log.success("✅ Using mkcert SSL certificates (trusted on this machine)");
      return { key: fs.readFileSync(keyPath), cert: fs.readFileSync(certPath) };
    }

    // No certificates found
    log.error("❌ No SSL certificates found in:", certDir);
    log.warning("💡 To generate trusted certificates, run:");
    log.warning("   cd backend/api/certs");
    log.warning(`   mkcert -key-file key.pem -cert-file cert.pem localhost 127.0.0.1 ::1 ${this.localIP}`);
    
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
      
      // Register with cleanup system if available
      if (cleanupRegistry) {
        cleanupRegistry.register(server);
      }
      
      return new Promise((resolve, reject) => {
        server.listen(this.actualPort, "0.0.0.0", () => {
          console.log(`🔒 HTTPS server running:
   Local: https://localhost:${this.actualPort}
   Network: https://${this.localIP}:${this.actualPort} (for mobile/other devices)
   ✅ No certificate warnings - trusted by mkcert`);
          resolve(server);
        });

        server.on("error", (error) => {
          if (error.code === "EADDRINUSE") {
            log.error(`❌ HTTPS port ${this.actualPort} is already in use`);
            log.warning("💡 HTTPS server will retry with a different port...");
            
            // Try to find another port and restart
            this.findAvailablePort(this.actualPort + 1)
              .then(newPort => {
                this.actualPort = newPort;
                log.warning(`🔄 Retrying HTTPS on port ${newPort}...`);
                server.close();
                this.start().then(resolve).catch(reject);
              })
              .catch(err => {
                log.error("❌ Could not find alternative HTTPS port:", err.message);
                resolve(null);
              });
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
