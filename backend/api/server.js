// server.js - Clean modular architecture
process.env.NODE_ENV ||= "development";

// Simple logger - only log on errors unless in debug mode
const isDebugMode = process.env.NODE_ENV === 'debug';
const log = {
  info: (msg) => isDebugMode && console.log(msg),
  success: (msg) => isDebugMode && console.log(msg),
  warning: (msg) => console.log(msg), // Keep warnings visible
  error: (msg) => console.error(msg)  // Always show errors
};

// ==================== IMPORTS ====================
const express = require("express");
const cors = require("cors");
const fs = require("fs");
const path = require("path");
const HttpsServer = require("./https-server");
const AtomPredictor = require("../services/AtomPredictor");
const MolecularProcessor = require("../services/molecular-processor");
// UserService import - only if database is available
let UserService = null;
try {
  UserService = require("../services/user-service");
} catch (error) {
  log.warning('⚠️ UserService not available - running without user management');
}
const {
  ImageMoleculeSchema,
  TextMoleculeSchema,
  SdfGenerationSchema,
} = require("../schemas/schemas");

// ==================== DATABASE CONFIGURATION ====================
let pool = null;
let dbConnected = false;

try {
  const { Pool } = require('pg');
  
  // Database configuration with local development defaults
  const dbConfig = {
    host: process.env.DB_HOST || 'localhost',
    port: process.env.DB_PORT || 5432,
    database: process.env.DB_NAME || 'mol_users',
    user: process.env.DB_USER || 'mol_user',
    password: process.env.DB_PASSWORD || 'mol_password',
    // Connection pool settings
    max: 20, // maximum number of clients in pool
    idleTimeoutMillis: 30000, // close idle clients after 30 seconds
    connectionTimeoutMillis: 2000, // return error after 2 seconds if connection could not be established
  };

  pool = new Pool(dbConfig);
  
  // Database connection error handling
  pool.on('error', (err, client) => {
    log.error('🔴 Unexpected error on idle client', err);
    log.info('💡 Database connection will be retried automatically');
  });
  
  log.success('✅ PostgreSQL module loaded successfully');
} catch (error) {
  log.warning('⚠️ PostgreSQL module not available - running without database');
  log.info('💡 Install with: npm install pg pg-pool');
}

// Database connection error handling (only if pool exists)
if (pool) {
  pool.on('error', (err, client) => {
    log.error('🔴 Unexpected error on idle client', err);
    log.info('💡 Database connection will be retried automatically');
  });
}

// Test database connection on startup
const testDatabaseConnection = async () => {
  if (!pool) {
    log.warning('⚠️ Database not available - running without persistent user storage');
    return false;
  }
  
  try {
    const client = await pool.connect();
    const result = await client.query('SELECT NOW()');
    client.release();
    log.success('✅ Database connected successfully');
    dbConnected = true;
    return true;
  } catch (err) {
    log.error('🔴 Database connection failed:', err.message);
    log.info('💡 Make sure PostgreSQL is running and credentials are correct');
    if (process.env.NODE_ENV === 'development') {
      log.info('💡 For local development, run: createdb mol_users');
    }
    dbConnected = false;
    return false;
  }
};

// ==================== CONFIGURATION ====================
const app = express();
const DEFAULT_PORT = 8080;
const HTTPS_PORT = process.env.HTTPS_PORT || 3001;

// Utility function to get local IP address
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

// Utility to attempt cleanup of processes using our ports
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

const findAvailablePort = async (startPort) => {
  const net = require("net");

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
    // Try up to 100 ports
    if (await isPortAvailable(port)) {
      return port;
    }
    port++;
  }
  throw new Error(
    `No available ports found between ${startPort} and ${startPort + 100}`,
  );
};

const PORT = process.env.PORT || DEFAULT_PORT;

// Initialize modules
const atomPredictor = new AtomPredictor(process.env.OPENAI_API_KEY);
const molecularProcessor = new MolecularProcessor();
const userService = (pool && UserService) ? new UserService(pool) : null;

// Initialize database on startup
const initializeDatabase = async () => {
  if (!userService) {
    console.log('⚠️ User service not available - running without persistent user storage');
    return;
  }
  
  try {
    const dbConnected = await testDatabaseConnection();
    if (dbConnected) {
      await userService.initializeTables();
      console.log('✅ Database initialized successfully');
    } else {
      console.log('⚠️ Database not connected - running without persistent user storage');
    }
  } catch (error) {
    console.error('🔴 Database initialization failed:', error.message);
    console.log('💡 Server will continue but user data will not persist');
  }
};

// ==================== MIDDLEWARE ====================
app.use(cors());
app.use(express.json({ limit: "50mb" }));

// ==================== HEALTH CHECK ENDPOINT ====================
app.get('/health', (req, res) => {
  res.json({ 
    status: 'ok', 
    timestamp: new Date().toISOString(),
    version: process.env.npm_package_version || '1.0.0'
  });
});

// ==================== ERROR LOGGING ENDPOINT ====================
app.post('/api/log-error', (req, res) => {
  const error = req.body;
  
  // Log to console for AI reading - STRUCTURED FORMAT
  console.error('🚨 FRONTEND ERROR:', JSON.stringify({
    timestamp: error.timestamp,
    type: error.type,
    message: error.message,
    source: error.source,
    userAgent: req.get('User-Agent'),
    ip: req.ip || req.connection.remoteAddress
  }, null, 2));
  
  // Optionally log to file or database here
  
  res.status(200).json({ status: 'logged' });
});

// User data now stored in PostgreSQL instead of in-memory
// Database schema will be created by the database setup script

// ==================== DEVELOPMENT MIDDLEWARE ====================
// Live reload enabled for local development
if (process.env.NODE_ENV === "development" || process.env.NODE_ENV === undefined) {
  const livereload = require("livereload");
  const connectLivereload = require("connect-livereload");
  const net = require("net");

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

  const startLiveReload = async (port) => {
    try {
      const available = await isPortAvailable(port);
      if (!available) {
        console.log(`LiveReload port ${port} is already in use`);
        if (port === 35729) {
          console.log("Trying alternative port 35730...");
          return startLiveReload(35730);
        } else {
          console.log("Continuing without LiveReload...");
          return;
        }
      }

      const liveReloadServer = livereload.createServer({
        exts: ["html", "css", "js"],
        ignore: ["node_modules/**", "sdf_files/**", "*.log"],
        port: port,
      });

      // Add error handler for unexpected errors
      liveReloadServer.server.on("error", (err) => {
        console.error("LiveReload server error:", err.message);
      });

      // Only proceed if server starts successfully
      liveReloadServer.server.once("listening", () => {
        const frontendPath = path.join(__dirname, "..", "..", "frontend");
        liveReloadServer.watch(frontendPath);
        app.use(connectLivereload());

        liveReloadServer.server.once("connection", () => {
          setTimeout(() => {
            liveReloadServer.refresh("/");
          }, 100);
        });

        console.log(`🔄 LiveReload server started on port ${port}`);
        console.log(`👀 Watching frontend files: ${frontendPath}`);
      });
    } catch (err) {
      console.log("LiveReload server failed to start:", err.message);
      console.log("Continuing without LiveReload...");
    }
  };

  // Start LiveReload with fallback port
  startLiveReload(35729);
}

app.use(express.static(path.join(__dirname, "..", "..", "frontend", "core")));
app.use("/assets", express.static(path.join(__dirname, "..", "..", "frontend", "assets")));
app.use("/components", express.static(path.join(__dirname, "..", "..", "frontend", "components")));
app.use("/sdf_files", express.static(path.join(__dirname, "..", "..", "data", "sdf_files")));

// ==================== PAYMENT ROUTES ====================

// Stripe configuration endpoint
app.get("/stripe-config", (req, res) => {
  res.json({
    publishableKey: process.env.STRIPE_PUBLISHABLE_KEY || 'pk_test_demo_key_for_development'
  });
});

// Setup payment method endpoint
app.post("/setup-payment-method", async (req, res) => {
  try {
    const { payment_method, device_info, name } = req.body;
    
    if (!payment_method || !device_info) {
      return res.status(400).json({ error: "Payment method and device info required" });
    }

    // Generate device token
    const deviceToken = Buffer.from(`${device_info}-${Date.now()}-${Math.random()}`).toString('base64').replace(/[+/=]/g, '');
    
    // Create user in database
    const userData = {
      deviceToken,
      paymentMethodId: payment_method,
      deviceInfo: device_info,
      name: name || null
    };
    
    if (!userService) {
      return res.status(503).json({ error: "User service not available" });
    }
    
    const user = await userService.createUser(userData);
    
    // In production, you would:
    // 1. Create customer in Stripe
    // 2. Attach payment method to customer
    // 3. Handle 3D Secure if needed
    // For demo, we'll just return success
    
    res.json({
      success: true,
      device_token: deviceToken,
      requires_action: false // Set to true if 3D Secure needed
    });
    
  } catch (error) {
    console.error("Payment setup error:", error);
    if (error.message === 'Device token already exists') {
      res.status(409).json({ error: "Device already registered" });
    } else {
      res.status(500).json({ error: "Failed to setup payment method" });
    }
  }
});

// Update payment method endpoint
app.post("/update-payment-method", async (req, res) => {
  try {
    const { device_token, payment_method, name } = req.body;
    
    if (!device_token || !payment_method) {
      return res.status(400).json({ error: "Device token and payment method required" });
    }
    
    if (!userService) {
      return res.status(503).json({ error: "User service not available" });
    }
    
    const user = await userService.getUserByDeviceToken(device_token);
    if (!user) {
      return res.status(404).json({ error: "User not found" });
    }
    
    // Update user in database
    await userService.updateUser(device_token, {
      paymentMethodId: payment_method,
      name: name || user.name
    });
    
    // In production, you would:
    // 1. Update payment method in Stripe
    // 2. Handle 3D Secure if needed
    
    res.json({
      success: true,
      message: "Payment method updated successfully"
    });
    
  } catch (error) {
    console.error("Payment update error:", error);
    res.status(500).json({ error: "Failed to update payment method" });
  }
});

// Get payment methods endpoint
app.get("/get-payment-methods", async (req, res) => {
  try {
    const { device_token } = req.query;
    
    if (!device_token) {
      return res.status(400).json({ error: "Device token required" });
    }
    
    if (!userService) {
      return res.status(503).json({ error: "User service not available" });
    }
    
    const user = await userService.getUserByDeviceToken(device_token);
    if (!user) {
      return res.status(404).json({ error: "User not found" });
    }
    
    // In production, you would:
    // 1. Fetch payment methods from Stripe
    // 2. Return masked card details
    
    // For demo, return mock card data
    const cardInfo = {
      id: user.payment_method_id || 'pm_demo',
      last4: '4242',
      brand: 'visa',
      exp_month: 12,
      exp_year: 2025,
      is_default: true
    };
    
    res.json({
      payment_methods: [cardInfo],
      default_method: cardInfo.id
    });
    
  } catch (error) {
    console.error("Get payment methods error:", error);
    res.status(500).json({ error: "Failed to get payment methods" });
  }
});

// Delete payment method endpoint
app.delete("/delete-payment-method", async (req, res) => {
  try {
    const { device_token, payment_method_id } = req.body;
    
    if (!device_token || !payment_method_id) {
      return res.status(400).json({ error: "Device token and payment method ID required" });
    }
    
    if (!userService) {
      return res.status(503).json({ error: "User service not available" });
    }
    
    const user = await userService.getUserByDeviceToken(device_token);
    if (!user) {
      return res.status(404).json({ error: "User not found" });
    }
    
    // In production, you would:
    // 1. Detach payment method from Stripe customer
    // 2. Update database accordingly
    
    // For demo, just clear the payment method
    await userService.updateUser(device_token, {
      paymentMethodId: null
    });
    
    res.json({
      success: true,
      message: "Payment method deleted successfully"
    });
    
  } catch (error) {
    console.error("Delete payment method error:", error);
    res.status(500).json({ error: "Failed to delete payment method" });
  }
});

// Set default payment method endpoint
app.post("/set-default-payment-method", async (req, res) => {
  try {
    const { device_token, payment_method_id } = req.body;
    
    if (!device_token || !payment_method_id) {
      return res.status(400).json({ error: "Device token and payment method ID required" });
    }
    
    if (!userService) {
      return res.status(503).json({ error: "User service not available" });
    }
    
    const user = await userService.getUserByDeviceToken(device_token);
    if (!user) {
      return res.status(404).json({ error: "User not found" });
    }
    
    // In production, you would:
    // 1. Set default payment method in Stripe
    // 2. Update database accordingly
    
    await userService.updateUser(device_token, {
      paymentMethodId: payment_method_id
    });
    
    res.json({
      success: true,
      message: "Default payment method updated successfully"
    });
    
  } catch (error) {
    console.error("Set default payment method error:", error);
    res.status(500).json({ error: "Failed to set default payment method" });
  }
});

// Validate payment method endpoint
app.post("/validate-payment", async (req, res) => {
  try {
    const { device_token } = req.body;
    
    if (!device_token) {
      return res.status(400).json({ error: "Device token required" });
    }
    
    if (!userService) {
      return res.status(503).json({ error: "User service not available" });
    }
    
    const user = await userService.getUserByDeviceToken(device_token);

    if (!user) {
      return res.status(404).json({ error: "User not found" });
    }
    
    res.json({
      valid: true,
      user: {
        name: user.name,
        usage: user.usage,
        device_token: user.device_token
      }
    });
    
  } catch (error) {
    console.error("Payment validation error:", error);
    res.status(500).json({ error: "Failed to validate payment" });
  }
});

// Increment usage endpoint
app.post("/increment-usage", async (req, res) => {
  try {
    const { device_token } = req.body;
    
    if (!device_token) {
      return res.status(400).json({ error: "Device token required" });
    }
    
    if (!userService) {
      return res.status(503).json({ error: "User service not available" });
    }
    
    const usage = await userService.incrementUsage(device_token);
    
    res.json({
      usage: usage
    });
    
  } catch (error) {
    console.error("Usage increment error:", error);
    if (error.message === 'User not found') {
      res.status(404).json({ error: "User not found" });
    } else {
      res.status(500).json({ error: "Failed to increment usage" });
    }
  }
});

// ==================== ANALYSIS ROUTES ====================

// Image analysis route
app.post("/image-molecules", async (req, res) => {
  try {
    // Validate input schema
    const validation = ImageMoleculeSchema.safeParse(req.body);
    if (!validation.success) {
      return res.status(400).json({
        error: "Invalid input data",
        details: validation.error.issues,
      });
    }

    const {
      imageBase64,
      croppedImageBase64,
      x,
      y,
      cropMiddleX,
      cropMiddleY,
      cropSize,
    } = req.body;

    if (!imageBase64) {
      return res.status(400).json({ error: "No image data provided" });
    }

    const result = await atomPredictor.analyzeImage(
      imageBase64,
      croppedImageBase64,
      x,
      y,
      cropMiddleX,
      cropMiddleY,
      cropSize,
    );
    res.json(result);
  } catch (error) {
    console.error("Image analysis error:", error);
    res.status(500).json({ error: error.message });
  }
});

// Text analysis route (frontend compatibility endpoint)
app.post("/analyze-text", async (req, res) => {
  try {
    const { object } = req.body;
    
    if (!object || typeof object !== "string" || object.trim().length === 0) {
      return res.status(400).json({ error: "Invalid object name: must be a non-empty string" });
    }

    const result = await atomPredictor.analyzeText(object || "");
    res.json(result);
  } catch (error) {
    console.error("Text analysis error:", error);
    res.status(500).json({
      error: `Analysis failed: ${error.message}`,
    });
  }
});

// Text analysis route (legacy endpoint)
app.post("/object-molecules", async (req, res) => {
  try {
    // Validate input schema
    const validation = TextMoleculeSchema.safeParse(req.body);
    if (!validation.success) {
      return res.status(400).json({
        error: "Invalid input data",
        details: validation.error.issues,
      });
    }

    const { object } = req.body;

    if (!object) {
      return res.status(400).json({ error: "No object description provided" });
    }

    const result = await atomPredictor.analyzeText(object);
    res.json(result);
  } catch (error) {
    console.error("Text analysis error:", error);

    // Provide more specific error messages
    let errorMessage = error.message;
    let statusCode = 500;

    if (error.message.includes("network") || error.message.includes("fetch")) {
      errorMessage =
        "Network error: Unable to connect to AI service. Please check your internet connection.";
      statusCode = 503;
    } else if (
      error.message.includes("API key") ||
      error.message.includes("authentication")
    ) {
      errorMessage = "Authentication error: Invalid or missing API key.";
      statusCode = 401;
    } else if (
      error.message.includes("rate limit") ||
      error.message.includes("quota")
    ) {
      errorMessage = "Rate limit exceeded: Please try again later.";
      statusCode = 429;
    } else if (error.message.includes("timeout")) {
      errorMessage =
        "Request timeout: The AI service is taking too long to respond.";
      statusCode = 408;
    }

    res.status(statusCode).json({
      error: errorMessage,
      details: process.env.NODE_ENV === "development" ? error.stack : undefined,
    });
  }
});

// SDF generation route
app.post("/generate-sdfs", async (req, res) => {
  try {
    // Validate input schema
    const validation = SdfGenerationSchema.safeParse(req.body);
    if (!validation.success) {
      return res.status(400).json({
        error: "Invalid input data",
        details: validation.error.issues,
      });
    }

    const { smiles, overwrite = false } = req.body;

    if (!smiles || !Array.isArray(smiles)) {
      return res.status(400).json({ error: "smiles array is required" });
    }

    const result = await molecularProcessor.processSmiles(smiles, overwrite);

    res.json({
      sdfPaths: result.sdfPaths,
      errors: result.errors,
      skipped: result.skipped,
      message: `Generated ${result.sdfPaths.length} 3D structures from ${smiles.length} SMILES`,
    });
  } catch (error) {
    console.error("SDF generation error:", error);

    // Provide more specific error messages
    let errorMessage = error.message;
    let statusCode = 500;

    if (error.message.includes("network") || error.message.includes("fetch")) {
      errorMessage =
        "Network error: Unable to connect to molecular service. Please check your internet connection.";
      statusCode = 503;
    } else if (
      error.message.includes("file system") ||
      error.message.includes("permission")
    ) {
      errorMessage =
        "File system error: Unable to create SDF files. Please check directory permissions.";
      statusCode = 500;
    } else if (error.message.includes("timeout")) {
      errorMessage =
        "Request timeout: Molecular structure generation is taking too long.";
      statusCode = 408;
    }

    res.status(statusCode).json({
      error: errorMessage,
      details: process.env.NODE_ENV === "development" ? error.stack : undefined,
    });
  }
});

// Static routes
app.get("/", (req, res) => {
  res.sendFile(path.join(__dirname, "..", "..", "frontend", "core", "index.html"));
});

// SEO routes
app.get("/robots.txt", (req, res) => {
  res.type('text/plain');
  res.sendFile(path.join(__dirname, "..", "..", "frontend", "core", "robots.txt"));
});

app.get("/sitemap.xml", (req, res) => {
  res.type('application/xml');
  res.sendFile(path.join(__dirname, "..", "..", "frontend", "core", "sitemap.xml"));
});

// Catch-all route for SPA (prevents 404s on direct navigation)
app.get("*", (req, res, next) => {
  // Skip API routes, assets, and special files
  if (req.url.startsWith("/api/") || 
      req.url.startsWith("/assets/") || 
      req.url.startsWith("/components/") ||
      req.url.startsWith("/sdf_files/") ||
      req.url.includes(".")) {
    return next();
  }
  
  // Serve index.html for all other routes (SPA behavior)
  res.sendFile(path.join(__dirname, "..", "..", "frontend", "core", "index.html"));
});

// Request logging middleware
app.use((req, res, next) => {
  // Skip logging Chrome DevTools discovery requests
  if (!req.url.includes(".well-known/appspecific/com.chrome.devtools")) {
    log.info(`Incoming request: ${req.method} ${req.url}`);
  }

  // Handle Chrome DevTools discovery request
  if (req.url === "/.well-known/appspecific/com.chrome.devtools.json") {
    return res.status(404).json({ error: "Not found" });
  }

  next();
});

// ==================== SERVER STARTUP ====================
const isCloudFunction =
  process.env.FUNCTION_NAME ||
  process.env.FUNCTION_TARGET ||
  process.env.K_SERVICE ||
  process.env.GOOGLE_CLOUD_PROJECT ||
  process.env.GCP_PROJECT;
const isNetlify = process.env.NETLIFY;
const isTestMode =
  process.env.NODE_ENV === "test" || process.env.JEST_WORKER_ID;
const isServerless = isCloudFunction || isNetlify;

// Store server instances for cleanup
let httpServer;
let httpsServerInstance;

if (!isServerless && !isTestMode) {
  // Local development mode
  const startServer = async () => {
    try {
      let actualPort = PORT;

      // If default port is in use, try cleanup and then find an available port
      if (PORT === DEFAULT_PORT) {
        try {
          // First check if port is available
          if (!(await (async () => {
            const net = require("net");
            return new Promise((resolve) => {
              const server = net.createServer();
              server.listen(PORT, "0.0.0.0", () => {
                server.once("close", () => resolve(true));
                server.close();
              });
              server.on("error", () => resolve(false));
            });
          })())) {
            console.log(`⚠️ Port ${PORT} is in use, attempting cleanup...`);
            const cleanupSuccessful = await attemptPortCleanup(PORT);
            
            if (cleanupSuccessful) {
              console.log(`✅ Port ${PORT} cleanup completed, retrying...`);
              // Give a moment for the port to be fully released
              await new Promise(resolve => setTimeout(resolve, 1000));
            }
          }
          
          actualPort = await findAvailablePort(PORT);
          if (actualPort !== PORT) {
            console.log(
              `⚠️  Port ${PORT} is in use, using port ${actualPort} instead`,
            );
          }
        } catch (error) {
          console.error(`❌ Could not find available port: ${error.message}`);
          console.log(`💡 Try: pkill -f "node.*server.js" or use a different port`);
          process.exit(1);
        }
      }

      const localIP = getLocalIPAddress();
      httpServer = app.listen(actualPort, "0.0.0.0", () => {
        log.success(`✅ HTTP server running on http://localhost:${actualPort}`);
        log.info(`📱 Mobile access: http://${localIP}:${actualPort}`);
        
        // Always log server ready for tests
        if (process.env.NODE_ENV === 'test') {
          console.log('Server running on port', actualPort);
        }
      });
    } catch (error) {
      console.error(`❌ Failed to start server: ${error.message}`);
      process.exit(1);
    }
  };

  startServer()
    .then(() => {
      // Handle port binding errors
      httpServer.on("error", (error) => {
        if (error.code === "EADDRINUSE") {
          console.error(`❌ Port ${PORT} is already in use`);
          console.log(`💡 Solutions:`);
          console.log(
            `   1. Kill existing process: pkill -f "node.*server.js"`,
          );
          console.log(`   2. Use different port: PORT=8081 npm start`);
          console.log(`   3. Check what's using the port: lsof -i :${PORT}`);
          process.exit(1);
        } else if (error.code === "EACCES") {
          console.error(`❌ Permission denied: Cannot bind to port ${PORT}`);
          console.log(`💡 Try using a port > 1024 or run with sudo`);
          process.exit(1);
        } else {
          console.error("❌ Server error:", error.message);
          console.log(`💡 Check your network configuration and try again`);
          process.exit(1);
        }
      });

      // Initialize database after server starts
      initializeDatabase();
    })
    .catch((error) => {
      console.error(`❌ Failed to start server: ${error.message}`);
      process.exit(1);
    });

  // Start HTTPS server for development
  if (process.env.NODE_ENV !== "production") {
    const startHttpsServer = async () => {
      try {
        const httpsServer = new HttpsServer(app, HTTPS_PORT);
        httpsServerInstance = await httpsServer.start();
        
        if (httpsServerInstance) {
          log.success("✅ HTTPS server started successfully");
          
          // Handle HTTPS server errors after startup
          httpsServerInstance.on("error", (error) => {
            console.error("❌ HTTPS server error after startup:", error.message);
            console.log("💡 HTTPS server will continue running if possible");
          });
        } else {
          log.warning("⚠️ HTTPS server not started - continuing with HTTP only");
        }
      } catch (error) {
        console.error("❌ Failed to start HTTPS server:", error.message);
        console.log("💡 Continuing with HTTP server only");
      }
    };

    // Start HTTPS server after a short delay to avoid port conflicts
    setTimeout(startHttpsServer, 1000);
  }
} else {
  // Serverless mode
  if (isCloudFunction) {
    console.log(`Running in Cloud Functions mode`);
  } else if (isNetlify) {
    console.log(`Running in Netlify mode`);
  }

  if (process.env.NODE_ENV === "production" && !process.env.OPENAI_API_KEY) {
    console.error(
      "OPENAI_API_KEY environment variable is required for production",
    );
    process.exit(1);
  }
}

// Graceful shutdown handling for nodemon restarts
if (!isServerless && !isTestMode) {
  const gracefulShutdown = (signal) => {
    console.log(`\n${signal} received: closing HTTP/HTTPS servers gracefully`);

    const closeServer = (server, name) => {
      return new Promise((resolve) => {
        if (server) {
          server.close(() => {
            console.log(`${name} server closed`);
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
      console.log("Shutdown complete");
      process.exit(0);
    });

    // Force exit after 5 seconds
    setTimeout(() => {
      console.error("Forced shutdown after timeout");
      process.exit(1);
    }, 5000);
  };

  // Listen for termination signals
  process.on("SIGTERM", () => gracefulShutdown("SIGTERM"));
  process.on("SIGINT", () => gracefulShutdown("SIGINT"));
  process.on("SIGUSR2", () => gracefulShutdown("SIGUSR2")); // Nodemon uses this
}

// Global error handlers
process.on("unhandledRejection", (reason, promise) => {
  console.error("❌ Unhandled Rejection at:", promise, "reason:", reason);
  console.log(
    "💡 This usually indicates a network or API error. Check your internet connection and API keys.",
  );
});

process.on("uncaughtException", (error) => {
  console.error("❌ Uncaught Exception:", error.message);
  console.log("💡 Application crashed. Check the error details above.");
  process.exit(1);
});

// Always export the app for Cloud Functions
module.exports = app;
