const https = require("https");
const fs = require("fs");
const path = require("path");
const os = require("os");
const { execSync } = require("child_process");
const net = require("net");

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

  generateSelfSignedCert() {
    const certDir = path.join(__dirname, "certs");
    const keyPath = path.join(certDir, "key.pem");
    const certPath = path.join(certDir, "cert.pem");
    const configPath = path.join(certDir, "openssl.conf");

    if (fs.existsSync(keyPath) && fs.existsSync(certPath)) {

      return { key: fs.readFileSync(keyPath), cert: fs.readFileSync(certPath) };
    }

    if (!fs.existsSync(certDir)) fs.mkdirSync(certDir, { recursive: true });



    try {
      const configContent = `[req]
distinguished_name = req_distinguished_name
req_extensions = v3_req
prompt = no

[req_distinguished_name]
C = US
ST = Development
L = Local
O = MolecularAnalysisApp
OU = Development
CN = localhost

[v3_req]
basicConstraints = CA:FALSE
keyUsage = digitalSignature, keyEncipherment, dataEncipherment
extendedKeyUsage = serverAuth, clientAuth
subjectAltName = @alt_names

[alt_names]
DNS.1 = localhost
DNS.2 = *.localhost
DNS.3 = mol.local
DNS.4 = *.mol.local
IP.1 = 127.0.0.1
IP.2 = 0.0.0.0
IP.3 = ${this.localIP}
IP.4 = ::1
`;

      fs.writeFileSync(configPath, configContent);

      execSync(`openssl genrsa -out ${keyPath} 2048`, { stdio: "pipe" });
      execSync(
        `openssl req -new -x509 -key ${keyPath} -out ${certPath} -days 3650 -config ${configPath}`,
        { stdio: "pipe" },
      );

      fs.unlinkSync(configPath);

      return { key: fs.readFileSync(keyPath), cert: fs.readFileSync(certPath) };
    } catch (err) {
      return null;
    }
  }

  async start() {
    try {
      const credentials = this.generateSelfSignedCert();
      if (!credentials) {
        return null;
      }

      // Find an available port
      try {
        this.actualPort = await this.findAvailablePort(this.requestedPort);
        
        if (this.actualPort !== this.requestedPort) {
        }
      } catch (error) {
        return null;
      }

      const server = https.createServer(credentials, this.app);
      
      return new Promise((resolve, reject) => {
        server.listen(this.actualPort, "0.0.0.0", () => {
          resolve(server);
        });

        server.on("error", (error) => {
          if (error.code === "EADDRINUSE") {
            
            // Try to find another port and restart
            this.findAvailablePort(this.actualPort + 1)
              .then(newPort => {
                this.actualPort = newPort;
                server.close();
                this.start().then(resolve).catch(reject);
              })
              .catch(err => {
                resolve(null);
              });
          } else if (error.code === "EACCES") {
            resolve(null);
          } else {
            resolve(null);
          }
        });
      });

    } catch (error) {
      return null;
    }
  }

  // Graceful shutdown
  async stop(server) {
    if (server) {
      return new Promise((resolve) => {
        server.close(() => {
          resolve();
        });
      });
    }
  }
}

module.exports = HttpsServer;
