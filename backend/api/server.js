const config = require('../config/env');
config.validateConfig();
const isDebugMode = config.NODE_ENV === 'debug';
const log = isDebugMode
  ? {
      info: (msg) => console.log(msg),
      success: (msg) => console.log(msg),
      warning: (msg) => console.warn(msg),
      error: (msg) => console.error(msg),
    }
  : {
      info: () => {},
      success: () => {},
      warning: () => {},
      error: (msg) => console.error(msg),
    };


const express = require("express");
const cors = require("cors");
const fs = require("fs");
const path = require("path");
let fetchLib = null;
try { fetchLib = require('node-fetch'); } catch (_) { fetchLib = null; }
const HttpsServer = require("./https-server");
const Structuralizer = require("../services/Structuralizer");
const MolecularProcessor = require("../services/molecular-processor");
const { resolveName, getPropertiesByCID } = require("../services/name-resolver");

let ScreenshotService = null;
if (config.NODE_ENV === 'development') {
  try {
    ScreenshotService = require("../services/screenshot-service");
  } catch (error) {
    log.warning('⚠️ Screenshot service not available in development mode');
  }
}


let proxy = null;
if (config.NODE_ENV === 'development') {
  try {
    proxy = require('http-proxy-middleware');
  } catch (error) {
    log.warning('⚠️ http-proxy-middleware not available - install with: npm install http-proxy-middleware');
  }
}

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

// Skip database setup entirely if disabled
if (!config.DB_ENABLED) {
  log.info('💡 Database disabled via DB_ENABLED=false - running without user storage');
} else {
  try {
    const { Pool } = require('pg');
  
  // Database configuration with local development defaults
  const dbConfig = {
    host: config.DB_HOST,
    port: config.DB_PORT,
    database: config.DB_NAME,
    user: config.DB_USER,
    password: config.DB_PASSWORD,
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
    // Database connection failed - continue without persistence
    dbConnected = false;
    return false;
  }
};


const app = express();
const DEFAULT_PORT = 8080;
// Consistent port configuration with cleanup
const HTTPS_PORT = process.env.HTTPS_PORT || 3001;
const HTTP_PORT = process.env.HTTP_PORT || 3002;

// Port cleanup utility
const cleanupPorts = async () => {
  const net = require("net");
  const portsToCheck = [HTTPS_PORT, HTTP_PORT, 3003, 3004, 3005];

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

const PORT = config.PORT;

// Initialize modules
// Do not inject a fake key; rely on real OPENAI_API_KEY or none
const openaiApiKey = config.OPENAI_API_KEY;
const structuralizer = new Structuralizer(openaiApiKey);
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

// ==================== CONFIGURATION ENDPOINT ====================
app.get('/api/config', (req, res) => {
  const paymentConfig = config.getPaymentConfig();
  
  res.json({
    payments: paymentConfig,
    environment: config.NODE_ENV,
    timestamp: new Date().toISOString()
  });
});

// ==================== ERROR LOGGING ENDPOINT ====================
app.post('/api/log-error', (req, res) => {
  const payload = req.body || {};
  const type = (payload.type || 'log').toLowerCase();
  const label = type === 'error' ? '❌ FRONTEND ERROR' : type === 'warn' ? '⚠️ FRONTEND WARN' : 'ℹ️ FRONTEND LOG';
  const logger = type === 'error' ? console.error : type === 'warn' ? console.warn : console.log;

  const timestamp = payload.timestamp || new Date().toISOString();
  const source = payload.source || '-';
  const userAgent = req.get('User-Agent') || '-';
  const ip = req.ip || req.connection.remoteAddress || '-';
  logger(`${label}: "${payload.message}" ts=${timestamp} src=${source} ua=${userAgent} ip=${ip}`);

  res.status(200).json({ status: 'logged' });
});

// User data now stored in PostgreSQL instead of in-memory
// Database schema will be created by the database setup script



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

    const result = await structuralizer.structuralize({
      object: "", // Let image detection determine the object
      imageBase64,
      x,
      y
    });
    res.json(result);
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
});

// Deprecated aliases (kept for compatibility): prefer /structuralize-*
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

    const result = await structuralizer.structuralizeText(object);
    res.json(result);
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
});

// Text analysis route (alias for backward compatibility)
app.post("/analyze-text", async (req, res) => {
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

    const result = await structuralizer.structuralizeText(object);
    res.json(result);
  } catch (error) {
    res.status(500).json({ error: error.message });
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
    res.status(500).json({ error: error.message });
  }
});

// ==================== DIRECT PUBCHEM SDF FETCHERS ====================
const fetchText = async (url) => {
  try {
    if (typeof fetch !== 'undefined') {
      const resp = await fetch(url);
      if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
      return await resp.text();
    }
  } catch (_) {}
  if (fetchLib) {
    const resp = await fetchLib(url);
    if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
    return await resp.text();
  }
  throw new Error('No fetch available');
};

const ensureSdfDir = () => {
  const dir = path.join(
    __dirname,
    "..",
    "..",
    process.env.NODE_ENV === "test" ? "test/sdf_files" : "backend/sdf_files"
  );
  if (!fs.existsSync(dir)) fs.mkdirSync(dir, { recursive: true });
  return dir;
};

const sanitizeName = (raw) => {
  try {
    const base = String(raw || '').trim();
    if (!base) return 'unknown';
    const cleaned = base.replace(/[^a-zA-Z0-9]+/g, '_').replace(/^_+|_+$/g, '');
    if (cleaned.length <= 64) return cleaned;
    const crypto = require('crypto');
    const hash = crypto.createHash('md5').update(base).digest('hex').slice(0, 12);
    return `${cleaned.slice(0, 32)}_${hash}`;
  } catch (_) {
    return 'unknown';
  }
};

const saveSdf = (filename, text) => {
  const dir = ensureSdfDir();
  const file = path.join(dir, filename);
  fs.writeFileSync(file, text, 'utf8');
  return `/sdf_files/${filename}`;
};

// Unified PubChem SDF fetch: select single highest-priority identifier (CID → name → SMILES)
const fetchPubchemSdf = async ({ cid, name, smiles, record_type = '3d' }) => {
  const rt = String(record_type).toLowerCase() === '2d' ? '2d' : '3d';
  let selected = null;
  if (cid != null) {
    selected = { type: 'cid', value: String(cid) };
  } else if (typeof name === 'string' && name.trim()) {
    selected = { type: 'name', value: name.trim() };
  } else if (typeof smiles === 'string' && smiles.trim()) {
    selected = { type: 'smiles', value: smiles.trim() };
  }
  if (!selected) {
    const err = new Error('no valid identifier provided');
    err.code = 'INVALID';
    throw err;
  }
  try {
    const enc = encodeURIComponent(selected.value);
    const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/${selected.type}/${enc}/SDF?record_type=${rt}`;
    const text = await fetchText(url);
    return { text, used: selected.type, rt, token: `${selected.type.toUpperCase()}_${selected.value}` };
  } catch (e) {
    const err = new Error(e?.message || 'identifier fetch failed');
    err.code = 'NOT_FOUND';
    throw err;
  }
};

app.post('/pubchem/sdf', async (req, res) => {
  try {
    const { cid, name, smiles, record_type = '3d' } = req.body || {};
    if (
      cid == null &&
      !(typeof name === 'string' && name.trim()) &&
      !(typeof smiles === 'string' && smiles.trim())
    ) {
      return res.status(400).json({ error: 'Provide cid, name, or smiles' });
    }
    const out = await fetchPubchemSdf({ cid, name, smiles, record_type });
    const safe = sanitizeName(`${out.token}_${out.rt}`);
    const sdfPath = saveSdf(`${safe}.sdf`, out.text);
    res.json({ sdfPath, status: 'ok', source: 'pubchem', used: out.used });
  } catch (error) {
    if (error && error.code === 'NOT_FOUND') {
      return res.status(404).json({ error: 'not found from provided identifiers' });
    }
    res.status(500).json({ error: error.message || 'fetch failed' });
  }
});


// Programmatic name → SMILES conversion (PubChem-backed)
app.post("/name-to-smiles", async (req, res) => {
  try {
    const { items } = req.body || {};
    if (!Array.isArray(items)) {
      return res.status(400).json({ error: "items array is required" });
    }
    const results = [];
    for (const it of items) {
      let name = it?.name || '';
      let cid = it?.cid ?? null;
      let smiles = null;
      try {
        if (cid) {
          const props = await getPropertiesByCID(cid);
          smiles = props?.smiles || null;
          if (!smiles) {
            const res = await resolveName(name);
            cid = res?.cid || cid;
            smiles = res?.smiles || null;
          }
          // prefer canonical name when available
          name = props?.title || props?.iupac || name;
        } else {
          const res = await resolveName(name);
          cid = res?.cid || null;
          smiles = res?.smiles || null;
          name = res?.title || res?.iupac || name;
        }
      } catch (_) {}
      results.push({ name, cid, smiles, status: smiles ? 'ok' : 'lookup_required' });
    }
    res.json({ molecules: results });
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
});

// Direct precise name → SDF (single name)
app.post("/name-to-sdf", async (req, res) => {
  try {
    const { name, overwrite = false } = req.body || {};
    if (typeof name !== 'string' || name.trim().length === 0) {
      return res.status(400).json({ error: "name is required" });
    }

    const byName = await molecularProcessor.generateSDFByName(name, overwrite);
    if (!byName || !byName.sdfPath) {
      return res.json({ name, sdfPath: null, status: 'lookup_required' });
    }

    res.json({ name: byName.name || name, sdfPath: byName.sdfPath, status: 'ok' });
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
});

// ==================== STATIC FILE SERVING ====================
app.use(express.static(path.join(__dirname, "..", "..", "frontend")));
app.use(express.static(path.join(__dirname, "..", "..", "public"))); // PWA manifest and icons
app.use("/dist", express.static(path.join(__dirname, "..", "..", "frontend", "dist")));
app.use("/assets", express.static(path.join(__dirname, "..", "..", "frontend", "assets")));
app.use("/components", express.static(path.join(__dirname, "..", "..", "frontend", "components")));
// Serve SDF files directory (uses test folder in test env, otherwise production folder)
app.use(
  "/sdf_files",
  express.static(
    path.join(
      __dirname,
      "..",
      "..",
      process.env.NODE_ENV === "test" ? "test/sdf_files" : "backend/sdf_files"
    )
  )
);




let screenshotService = null;
if (config.NODE_ENV === 'development' && ScreenshotService) {
  screenshotService = new ScreenshotService();
  log.info("🔧 Screenshot service enabled for development");
} else {
  log.info("📦 Screenshot service disabled for production");
}

// Development-only screenshot routes
if (config.NODE_ENV === 'development' && screenshotService) {
  // Serve screenshot files for LLM access
  app.use("/api/screenshot", express.static(path.join(__dirname, "..", "..")));

  // Capture current app state
  app.post("/api/capture-screenshot", async (req, res) => {
    try {
      const { filename } = req.body;
      const result = await screenshotService.captureApp(filename);
      
      res.json({
        success: true,
        screenshot: result,
        message: "Screenshot captured successfully"
      });
    } catch (error) {
      log.error(`Screenshot capture failed: ${error.message}`);
      res.status(500).json({
        success: false,
        error: error.message
      });
    }
  });

  // Capture screenshot with text input
  app.post("/api/capture-with-input", async (req, res) => {
    try {
      const { inputText, filename } = req.body;
      
      if (!inputText) {
        return res.status(400).json({
          success: false,
          error: "inputText is required"
        });
      }
      
      const result = await screenshotService.captureWithInput(inputText, filename);
      
      res.json({
        success: true,
        screenshot: result,
        message: "Screenshot with input captured successfully"
      });
    } catch (error) {
      log.error(`Screenshot with input failed: ${error.message}`);
      res.status(500).json({
        success: false,
        error: error.message
      });
    }
  });

  // Capture screenshot after analysis
  app.post("/api/capture-analysis", async (req, res) => {
    try {
      const { inputText, filename } = req.body;
      
      if (!inputText) {
        return res.status(400).json({
          success: false,
          error: "inputText is required"
        });
      }
      
      const result = await screenshotService.captureAnalysis(inputText, filename);
      
      res.json({
        success: true,
        screenshot: result,
        message: "Analysis screenshot captured successfully"
      });
    } catch (error) {
      log.error(`Analysis screenshot failed: ${error.message}`);
      res.status(500).json({
        success: false,
        error: error.message
      });
    }
  });

  // List all available screenshots
  app.get("/api/screenshots", async (req, res) => {
    try {
      const screenshots = await screenshotService.listScreenshots();
      
      res.json({
        success: true,
        screenshots: screenshots,
        count: screenshots.length
      });
    } catch (error) {
      log.error(`Failed to list screenshots: ${error.message}`);
      res.status(500).json({
        success: false,
        error: error.message
      });
    }
  });

  // Get specific screenshot info
  app.get("/api/screenshot-info/:filename", async (req, res) => {
    try {
      const { filename } = req.params;
      const screenshotPath = await screenshotService.getScreenshotPath(filename);
      
      const stats = fs.statSync(screenshotPath);
      
      res.json({
        success: true,
        screenshot: {
          filename: filename,
          path: screenshotPath,
          url: `/api/screenshot/${filename}`,
          size: stats.size,
          created: stats.birthtime,
          modified: stats.mtime
        }
      });
    } catch (error) {
      res.status(404).json({
        success: false,
        error: error.message
      });
    }
  });

  // Cleanup old screenshots
  app.post("/api/cleanup-screenshots", async (req, res) => {
    try {
      const { keepCount = 10 } = req.body;
      await screenshotService.cleanupOldScreenshots(keepCount);
      
      res.json({
        success: true,
        message: `Cleaned up old screenshots, kept ${keepCount} most recent`
      });
    } catch (error) {
      log.error(`Screenshot cleanup failed: ${error.message}`);
      res.status(500).json({
        success: false,
        error: error.message
      });
    }
  });

  log.info("🔧 Screenshot routes enabled for development");
} else {
  // Stub routes for production
  const screenshotDisabled = (req, res) => {
    res.status(503).json({ 
      error: "Screenshot functionality disabled in production" 
    });
  };
  
  app.post("/api/capture-screenshot", screenshotDisabled);
  app.post("/api/capture-with-input", screenshotDisabled);
  app.post("/api/capture-analysis", screenshotDisabled);
  app.get("/api/screenshots", screenshotDisabled);
  app.get("/api/screenshot-info/:filename", screenshotDisabled);
  app.post("/api/cleanup-screenshots", screenshotDisabled);
  
  log.info("📦 Screenshot routes disabled for production");
}



// Stripe configuration endpoint
app.get("/api/stripe-config", (req, res) => {
  res.json({
    publishableKey: process.env.STRIPE_PUBLISHABLE_KEY || 'pk_test_demo_key_for_development'
  });
});

// Setup payment method endpoint
app.post("/api/setup-payment-method", async (req, res) => {
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
app.post("/api/update-payment-method", async (req, res) => {
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
app.get("/api/get-payment-methods", async (req, res) => {
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
app.post("/api/set-default-payment-method", async (req, res) => {
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

    const result = await structuralizer.structuralize({
      object: "", // Let image detection determine the object
      imageBase64,
      x,
      y
    });
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

    const result = await structuralizer.structuralizeText(object || "");
    res.json(result);
  } catch (error) {
    console.error("Text analysis error:", error);
    res.status(500).json({
      error: `Analysis failed: ${error.message}`,
    });
  }
});

// Preferred naming alias: structures-from-text
app.post("/structures-from-text", async (req, res) => {
  try {
    const { object } = req.body;
    if (!object || typeof object !== "string" || object.trim().length === 0) {
      return res.status(400).json({ error: "Invalid object name: must be a non-empty string" });
    }
    const result = await structuralizer.structuralizeText(object || "");
    res.json(result);
  } catch (error) {
    res.status(500).json({ error: `Analysis failed: ${error.message}` });
  }
});

// Structuralize alias routes (preferred naming)
app.post("/structuralize-text", async (req, res) => {
  try {
    const { object } = req.body;
    if (!object || typeof object !== "string" || object.trim().length === 0) {
      return res.status(400).json({ error: "Invalid object name: must be a non-empty string" });
    }
    const result = await structuralizer.structuralizeText(object || "");
    res.json(result);
  } catch (error) {
    res.status(500).json({ error: `Structuralization failed: ${error.message}` });
  }
});

app.post("/structuralize-image", async (req, res) => {
  try {
    const { imageBase64, croppedImageBase64, x, y, cropMiddleX, cropMiddleY, cropSize } = req.body;
    if (!imageBase64) {
      return res.status(400).json({ error: "No image data provided" });
    }
    const result = await structuralizer.structuralizeImage(
      imageBase64,
      croppedImageBase64,
      typeof x === 'number' ? x : null,
      typeof y === 'number' ? y : null,
      typeof cropMiddleX === 'number' ? cropMiddleX : null,
      typeof cropMiddleY === 'number' ? cropMiddleY : null,
      typeof cropSize === 'number' ? cropSize : null
    );
    res.json(result);
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
});

// Unified multimodal structuralization: accepts { object?, imageBase64? }
app.post("/structuralize", async (req, res) => {
  try {
    if (!req.body || (typeof req.body !== 'object')) {
      return res.status(400).json({ error: "Invalid payload" });
    }
    const out = await structuralizer. structuralize(req.body);
    res.json(out);
  } catch (error) {
    res.status(500).json({ error: `Structuralization failed: ${error.message}` });
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

    const result = await structuralizer.structuralizeText(object);
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

console.log('NODE_ENV:', config.NODE_ENV);
console.log('Setting up frontend routes');

// Helper to serve index.html with token replacement
const serveIndexWithTokens = (req, res) => {
  try {
    const indexPath = path.join(__dirname, "..", "..", "frontend", "core", "index.html");
    let html = fs.readFileSync(indexPath, "utf8");

    const host = req.get('host') || 'localhost';
    const baseUrl = `${req.protocol || 'http'}://${host}`;

    // Generate cache busting timestamp
    const cacheBust = Date.now();

    const tokens = {
      '{{TITLE}}': 'Molecular Contents',
      '{{DESCRIPTION}}': 'Advanced 3D visualization tool for chemical compounds. Analyze molecular structures from camera images or text input.',
      '{{CANONICAL}}': `${baseUrl}/`,
      '{{OG_IMAGE}}': '/assets/favicon.svg',
      '{{CACHE_BUST}}': cacheBust
    };

    html = html.replace(/\{\{TITLE\}\}|\{\{DESCRIPTION\}\}|\{\{CANONICAL\}\}|\{\{OG_IMAGE\}\}|\{\{CACHE_BUST\}\}/g, (match) => tokens[match] || '');

    res.type('html').send(html);
  } catch (e) {
    res.sendFile(path.join(__dirname, "..", "..", "frontend", "core", "index.html"));
  }
};

// Serve manifest.json with cache busting
app.get("/manifest.json", (req, res) => {
  try {
    const manifestPath = path.join(__dirname, "..", "..", "public", "manifest.json");
    let manifest = fs.readFileSync(manifestPath, "utf8");

    // Generate cache busting timestamp
    const cacheBust = Date.now();

    // Replace cache busting tokens
    manifest = manifest.replace(/\{\{CACHE_BUST\}\}/g, cacheBust);

    res.type('json').send(manifest);
  } catch (e) {
    res.status(500).send('Error serving manifest');
  }
});

app.get("/", serveIndexWithTokens);


app.get("/robots.txt", (req, res) => {
  res.type('text/plain');
  res.send(`User-agent: *
Allow: /

Sitemap: https://${req.get('host') || 'localhost'}/sitemap.xml
`);
});

app.get("/sitemap.xml", (req, res) => {
  res.type('application/xml');
  res.send(`<?xml version="1.0" encoding="UTF-8"?>
<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">
  <url>
    <loc>https://${req.get('host') || 'localhost'}/</loc>
    <changefreq>weekly</changefreq>
    <priority>1.0</priority>
  </url>
</urlset>
`);
});


app.use(express.static(path.join(__dirname, "..", "..", "frontend", "dist")));

// Catch-all route for SPA (prevents 404s on direct navigation)
app.get("*", (req, res, next) => {
  // Skip API routes, assets, and special files
  if (req.url.startsWith("/api/") || 
      req.url.startsWith("/assets/") || 
      req.url.startsWith("/dist/") ||
      req.url.startsWith("/components/") ||
      req.url.startsWith("/sdf_files/") ||
      req.url.includes(".")) {
    return next();
  }
  
  // Serve index.html for all other routes (SPA behavior)
  serveIndexWithTokens(req, res);
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
  config.FUNCTION_NAME ||
  config.FUNCTION_TARGET ||
  config.K_SERVICE ||
  config.GOOGLE_CLOUD_PROJECT ||
  config.GCP_PROJECT;
const isNetlify = !!config.NETLIFY;
const isTestMode =
  config.NODE_ENV === "test" || !!process.env.JEST_WORKER_ID;
const isIntegrationTest = !!config.INTEGRATION_TEST;
const isServerless = !!(isCloudFunction || isNetlify);

// Store server instances for cleanup
let httpServer;
let httpsServerInstance;

console.log('Server conditions:', { 
  isServerless: isServerless || false, 
  isTestMode: isTestMode || false, 
  isIntegrationTest: isIntegrationTest || false 
});

if (!isServerless && (!isTestMode || isIntegrationTest)) {
  console.log('Starting server in local development mode...');
  // Local development mode
  const startServer = async () => {
    // Clean up any stale ports before starting
    console.log('🧹 Cleaning up ports...');
    await cleanupPorts();

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
        console.log(`http://localhost:${actualPort}`);
        if (localIP) {
          console.log(`On your phone: http://${localIP}:${actualPort}`);
        }
      });
    } catch (error) {
      console.error(`❌ Failed to start server: ${error.message}`);
      process.exit(1);
    }
  };

  // Start HTTPS server with mkcert certificates
  const startHttpsServer = async () => {
      try {
        const httpsServer = new HttpsServer(app, HTTPS_PORT);
        httpsServerInstance = await httpsServer.start();
        
        if (httpsServerInstance) {
          log.success("✅ HTTPS server started with trusted certificates");

          // Handle HTTPS server errors after startup
          httpsServerInstance.on("error", (error) => {
            console.error("❌ HTTPS server error after startup:", error.message);
            console.log("💡 HTTPS server will continue running if possible");
          });

          // Register with cleanup system if available
          try {
            const cleanupRegistry = require('../test/fixtures/cleanup-registry');
            cleanupRegistry.register(httpsServerInstance);
          } catch (e) {
            // Cleanup registry not available
          }
        } else {
          log.warning("⚠️ HTTPS server not started - continuing with HTTP only");
        }
      } catch (error) {
        console.error("❌ Failed to start HTTPS server:", error.message);
        console.log("💡 Continuing with HTTP server only");
      }
    };

  // Prefer HTTPS: try it first; fall back to HTTP only if HTTPS not available
  startHttpsServer()
    .then(() => {
      if (!httpsServerInstance) {
        return startServer()
          .then(() => {
            // Handle port binding errors
            httpServer.on("error", (error) => {
              if (error.code === "EADDRINUSE") {
                console.error(`❌ Port ${PORT} is already in use`);
                console.log(`💡 Solutions:`);
                console.log(`   1. Kill existing process: pkill -f "node.*server.js"`);
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
      } else {
        // HTTPS started successfully; initialize database
        initializeDatabase();
      }
    })
    .catch((error) => {
      console.error(`❌ Failed to start HTTPS/HTTP server: ${error.message}`);
      process.exit(1);
    });
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

    // Force exit after 5 seconds
    setTimeout(() => {
      process.exit(1);
    }, 5000);
  };

  // Listen for termination signals
  process.on("SIGTERM", () => gracefulShutdown("SIGTERM"));
  process.on("SIGINT", () => gracefulShutdown("SIGINT"));
  process.on("SIGUSR2", () => gracefulShutdown("SIGUSR2")); // Nodemon uses this
}

// Global error handlers
const ErrorStackLogger = require("../services/ErrorStackLogger");
process.on("unhandledRejection", (reason) => {
  ErrorStackLogger.log(reason, "Unhandled Rejection");
});

process.on("uncaughtException", (error) => {
  ErrorStackLogger.log(error, "Uncaught Exception");
  process.exit(1);
});

// Hot reload is now handled by LiveReload in development

// Always export the app for Cloud Functions
module.exports = app;
// Also provide a named export for platforms expecting a handler name
module.exports.main = app;
module.exports.molecularAnalysis = app;
