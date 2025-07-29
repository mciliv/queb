/**
 * MOLECULAR ANALYSIS API SERVER
 * Purpose: Express.js API server for molecular analysis and visualization
 * Features: AI-powered analysis, 3D molecule generation, user management, payment processing
 * 
 * Service Architecture:
 * - AtomPredictor: OpenAI-based molecular analysis
 * - MolecularProcessor: SMILES to 3D structure conversion
 * - UserService: Account management (PostgreSQL or in-memory fallback)
 * - ErrorHandler: Centralized error handling and logging
 */

process.env.NODE_ENV ||= "development";

// ==================== CORE DEPENDENCIES ====================
const express = require("express");
const cors = require("cors");
const fs = require("fs");
const path = require("path");

// ==================== SERVICE IMPORTS ====================
const HttpsServer = require("./https-server");
const AtomPredictor = require("../services/AtomPredictor");
const MolecularProcessor = require("../services/molecular-processor");
const ErrorHandler = require("../services/error-handler");

// ==================== USER SERVICE INITIALIZATION ====================
// Service Priority: Full PostgreSQL UserService > Simple in-memory > No user management
let UserService = null;
let SimpleUserService = null;
let availableUserService = null;

try {
  // Attempt to load full PostgreSQL user service
  UserService = require("../services/user-service");
  console.log('✅ PostgreSQL UserService available');
} catch (error) {
  console.log('⚠️ PostgreSQL UserService not available, trying in-memory fallback');
  
  try {
    // Fallback to simple in-memory user service  
    SimpleUserService = require("../services/simple-user-service");
    console.log('✅ In-memory UserService loaded (no data persistence)');
  } catch (fallbackError) {
    console.log('⚠️ No user service available - server will run without user management');
    console.log('💡 User features will be disabled but molecular analysis will work');
  }
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
    console.error('🔴 Unexpected error on idle client', err);
    console.log('💡 Database connection will be retried automatically');
  });
  
  console.log('✅ PostgreSQL module loaded successfully');
} catch (error) {
  console.log('⚠️ PostgreSQL module not available - running without database');
  console.log('💡 Install with: npm install pg pg-pool');
}

// Database connection error handling (only if pool exists)
if (pool) {
  pool.on('error', (err, client) => {
    console.error('🔴 Unexpected error on idle client', err);
    console.log('💡 Database connection will be retried automatically');
  });
}

// Test database connection on startup
const testDatabaseConnection = async () => {
  if (!pool) {
    console.log('⚠️ Database not available - running without persistent user storage');
    return false;
  }
  
  try {
    const client = await pool.connect();
    const result = await client.query('SELECT NOW()');
    client.release();
    console.log('✅ Database connected successfully');
    dbConnected = true;
    return true;
  } catch (err) {
    console.error('🔴 Database connection failed:', err.message);
    console.log('💡 Make sure PostgreSQL is running and credentials are correct');
    if (process.env.NODE_ENV === 'development') {
      console.log('💡 For local development, run: createdb mol_users');
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



const PORT = process.env.PORT || DEFAULT_PORT;

// Initialize modules
const atomPredictor = new AtomPredictor(process.env.OPENAI_API_KEY);
const molecularProcessor = new MolecularProcessor();
// Initialize user service with detailed logging
let userService = null;
if (pool && UserService) {
  userService = new UserService(pool);
  console.log('✅ PostgreSQL UserService initialized');
} else if (SimpleUserService) {
  userService = new SimpleUserService();
  console.log('✅ Simple UserService initialized');
} else {
  console.log('⚠️ No user service available');
}

// Initialize database on startup
const initializeDatabase = async () => {
  if (!userService) {
    console.log('⚠️ User service not available - running without user storage');
    return;
  }
  
  try {
    // If using PostgreSQL service, test connection first
    if (pool && UserService && userService instanceof UserService) {
      const dbConnected = await testDatabaseConnection();
      if (dbConnected) {
        await userService.initializeTables();
        console.log('✅ PostgreSQL database initialized successfully');
      } else {
        console.log('⚠️ Database not connected - running without persistent user storage');
      }
    } else {
      // Using simple in-memory service
      await userService.initializeTables();
      console.log('✅ In-memory user service initialized successfully');
    }
  } catch (error) {
    console.error('🔴 User service initialization failed:', error.message);
    console.log('💡 Server will continue but user data may not persist');
  }
};

// ==================== MIDDLEWARE ====================
app.use(cors());
app.use(express.json({ limit: "50mb" }));

// LiveReload for development
if (process.env.NODE_ENV === 'development') {
  try {
    const livereload = require('livereload');
    const connectLivereload = require('connect-livereload');
    
    // Create livereload server
    const liveReloadServer = livereload.createServer({
      exts: ['html', 'css', 'js'],
      port: 35730
    });
    
    // Watch frontend directory
    liveReloadServer.watch(path.join(__dirname, '../../frontend'));
    
    // Add livereload middleware
    app.use(connectLivereload({ port: 35730 }));
    
    console.log('🔄 LiveReload server started on port 35730');
  } catch (error) {
    console.log('⚠️ LiveReload not available:', error.message);
  }
}

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

// Logging endpoint for frontend logs
app.post('/api/logs', (req, res) => {
  try {
    const { logs, timestamp, userAgent, url } = req.body;
    
    if (!logs || !Array.isArray(logs)) {
      return res.status(400).json({ error: 'Invalid logs format' });
    }

    // Process each log entry
    logs.forEach(logEntry => {
      const { level, message, data, timestamp: logTimestamp, source } = logEntry;
      
      // Format the log message
      const timestampStr = new Date(logTimestamp).toISOString();
      const sourceStr = source ? `[${source.file}:${source.line}]` : '';
      const levelStr = level.toUpperCase().padEnd(5);
      
      let logMessage = `[${timestampStr}] ${levelStr} ${sourceStr} ${message}`;
      
      // Add data if present
      if (data) {
        if (typeof data === 'object') {
          logMessage += ` | ${JSON.stringify(data)}`;
        } else {
          logMessage += ` | ${data}`;
        }
      }
      
      // Output to appropriate stream based on level
      switch (level) {
        case 'error':
          console.error(logMessage);
          break;
        case 'warn':
          console.warn(logMessage);
          break;
        case 'debug':
          console.debug(logMessage);
          break;
        default:
          console.log(logMessage);
      }
    });

    // Log request metadata
    console.log(`[${new Date().toISOString()}] INFO [LOGGER] Received ${logs.length} logs from ${url} (${userAgent})`);

    res.json({ success: true, processed: logs.length });
  } catch (error) {
    console.error('Error processing logs:', error);
    res.status(500).json({ error: 'Failed to process logs' });
  }
});

// User data now stored in PostgreSQL instead of in-memory
// Database schema will be created by the database setup script

// ==================== STATIC FILES ====================

// Disable caching for HTML files in development
if (process.env.NODE_ENV === "development") {
  app.use(express.static(path.join(__dirname, "..", "..", "frontend", "core"), {
    setHeaders: (res, path) => {
      if (path.endsWith('.html')) {
        res.set('Cache-Control', 'no-cache, no-store, must-revalidate');
        res.set('Pragma', 'no-cache');
        res.set('Expires', '0');
      }
    }
  }));
} else {
  app.use(express.static(path.join(__dirname, "..", "..", "frontend", "core")));
}
app.use("/assets", express.static(path.join(__dirname, "..", "..", "frontend", "assets")));

// Favicon route - serve SVG as favicon.ico
app.get('/favicon.ico', (req, res) => {
  if (process.env.NODE_ENV === "development") {
    res.set('Cache-Control', 'no-cache, no-store, must-revalidate');
    res.set('Pragma', 'no-cache');
    res.set('Expires', '0');
  }
  res.type('image/svg+xml');
  res.sendFile(path.join(__dirname, "..", "..", "frontend", "assets", "favicon.svg"));
});

// Disable caching for components in development to prevent issues with cached JavaScript
if (process.env.NODE_ENV === "development") {
  app.use("/components", express.static(path.join(__dirname, "..", "..", "frontend", "components"), {
    setHeaders: (res) => {
      res.set('Cache-Control', 'no-cache, no-store, must-revalidate');
      res.set('Pragma', 'no-cache');
      res.set('Expires', '0');
    }
  }));
} else {
  app.use("/components", express.static(path.join(__dirname, "..", "..", "frontend", "components")));
}

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
      console.log('🔴 UserService is null during payment setup');
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
    res.json({ output: result });
  } catch (error) {
    console.error("Image analysis error:", error);
    res.status(500).json({ error: error.message });
  }
});

// Text analysis route (frontend compatibility endpoint)
app.post("/analyze-text", async (req, res) => {
  try {
    const { text } = req.body;
    
    if (!text) {
      return res.status(400).json({ error: "No text provided" });
    }

    // Map to the existing object-molecules logic
    const result = await atomPredictor.analyzeText(text);
    res.json({ output: result });
  } catch (error) {
    const { errorMessage, statusCode } = ErrorHandler.handleAIError(error, 'text analysis endpoint');
    res.status(statusCode).json({
      error: errorMessage,
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
    res.json({ output: result });
  } catch (error) {
    const { errorMessage, statusCode } = ErrorHandler.handleAIError(error, 'object-molecules endpoint');
    res.status(statusCode).json({
      error: errorMessage,
      details: process.env.NODE_ENV === "development" ? error.stack : undefined,
    });
  }
});


// SDF generation route (returns JSON with paths)
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
  res.sendFile(path.join(__dirname, "index.html"));
});

// Request logging middleware
app.use((req, res, next) => {
  // Skip logging Chrome DevTools discovery requests
  if (!req.url.includes(".well-known/appspecific/com.chrome.devtools")) {
    console.log(`Incoming request: ${req.method} ${req.url}`);
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
  // Simple server startup
  const localIP = getLocalIPAddress();
  
  httpServer = app.listen(PORT, "0.0.0.0", () => {
    console.log(`✅ HTTP server running on http://localhost:${PORT}`);
    console.log(`📱 Mobile access: http://${localIP}:${PORT}`);
  });

  // Handle port binding errors
  httpServer.on("error", (error) => {
    if (error.code === "EADDRINUSE") {
      console.error(`❌ Port ${PORT} is already in use`);
      console.log(`💡 Try: ./dev (kills existing processes)`);
      process.exit(1);
    } else {
      console.error("❌ Server error:", error.message);
      process.exit(1);
    }
  });

  // Initialize database after server starts
  initializeDatabase();

  // Start HTTPS server for development
  if (process.env.NODE_ENV !== "production") {
    setTimeout(async () => {
      try {
        const httpsServer = new HttpsServer(app, HTTPS_PORT);
        httpsServerInstance = await httpsServer.start();
        if (httpsServerInstance) {
          console.log("✅ HTTPS server started successfully");
        }
      } catch (error) {
        console.log("⚠️ HTTPS server not started - continuing with HTTP only");
      }
    }, 1000);
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
