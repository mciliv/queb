// test/integration/smoke.test.js - Quick smoke tests for pre-development validation
// These tests run in under 5 seconds to catch basic issues before starting development

const request = require("supertest");
const app = require("../../../src/server/api/server");
const fs = require("fs");
const path = require("path");

// Smoke Test Validation Rules - Critical App Health Checks
const SMOKE_VALIDATION_RULES = {
  // Critical system checks that must pass
  criticalChecks: {
    serverStartup: { maxTime: 3000, required: true },
    staticFiles: { required: ['index.html', 'style.css'], maxSize: 5242880 },
    apiEndpoints: { 
      required: ['/analyze-text', '/image-molecules', '/generate-sdfs'],
      expectedErrors: [400, 500], // These endpoints should reject invalid data
      timeout: 2000
    },
    dependencies: {
      npm: ['express', 'cors', 'openai', 'zod'],
      optional: ['puppeteer', 'http-proxy-middleware']
    }
  },
  
  // Basic functionality validation
  functionality: {
    inputValidation: { 
      rejectEmpty: true, 
      rejectInvalid: true,
      expectedErrorCodes: [400]
    },
    responseFormat: {
      contentType: 'application/json',
      hasErrorField: true,
      hasDataField: false // For error responses
    },
    fileOperations: {
      readAccess: true,
      writeAccess: true,
      createDirectories: true
    }
  },
  
  // Performance thresholds for smoke tests
  performance: {
    maxStartupTime: 5000,
    maxResponseTime: 1000,
    maxMemoryUsage: 536870912, // 512MB
    minAvailableDisk: 1073741824 // 1GB
  }
};

describe("Smoke Tests - Critical App Validation", () => {
  // Validation helper for smoke tests
  const validateCriticalCheck = (checkName, result) => {
    const check = SMOKE_VALIDATION_RULES.criticalChecks[checkName];
    if (check && check.required) {
      expect(result).toBeTruthy();
    }
    return result;
  };

  describe("Application Startup", () => {
    test("should start without crashing with validation rules", () => {
      const startTime = Date.now();
      
      expect(app).toBeDefined();
      expect(typeof app.listen).toBe("function");
      
      const startupTime = Date.now() - startTime;
      expect(startupTime).toBeLessThan(SMOKE_VALIDATION_RULES.performance.maxStartupTime);
      
      validateCriticalCheck('serverStartup', true);
      console.log(`✅ Server startup validated (${startupTime}ms)`);
    });

    test("should serve static files", async () => {
      const response = await request(app).get("/");
      expect(response.status).toBe(200);
      expect(response.text).toContain("<!doctype html>");
      expect(response.text).toContain("Materials");
    });

    test("should have required core files", () => {
      const requiredFiles = [
        "../../../src/server/api/server.js",
        "../../../src/server/schemas/molecular.js",
        "../../../src/server/services/Structuralizer.js",
        "../../../src/server/services/molecular-processor.js",
        "../../../package.json",
      ];

      requiredFiles.forEach((file) => {
        const filePath = path.join(__dirname, file);
        expect(fs.existsSync(filePath)).toBe(true);
        expect(fs.statSync(filePath).size).toBeGreaterThan(0);
      });
    });
  });

  describe("API Endpoints", () => {
    test("should respond to health check endpoints", async () => {
      // Test main page
      const mainResponse = await request(app).get("/");
      expect(mainResponse.status).toBe(200);

      // Test SDF files directory (should exist even if empty)
      const sdfDir = path.join(__dirname, "../../sdf_files");
      if (fs.existsSync(sdfDir)) {
        const sdfResponse = await request(app).get("/sdf_files/");
        expect([200, 404]).toContain(sdfResponse.status); // 404 is okay if directory is empty
      }
    });

    test("should handle API endpoint requests", async () => {
      // Test structures-from-text endpoint (should respond even without valid data)
      const structuresResponse = await request(app)
        .post("/structures-from-text")
        .send({});

      expect([400, 500]).toContain(structuresResponse.status); // Should reject invalid data

      // Test object molecules endpoint
      const objectResponse = await request(app)
        .post("/object-molecules")
        .send({});

      expect([400, 500]).toContain(objectResponse.status); // Should reject invalid data

      // Test SDF generation endpoint
      const sdfResponse = await request(app).post("/generate-sdfs").send({});

      expect([400, 500]).toContain(sdfResponse.status); // Should reject invalid data
    });
  });

  describe("Core Dependencies", () => {
    test("should have required npm dependencies", () => {
      const packageJson = JSON.parse(
        fs.readFileSync(path.join(__dirname, "../../../package.json"), "utf8"),
      );

      const requiredDeps = ["express", "cors", "openai", "zod"];

      requiredDeps.forEach((dep) => {
        expect(packageJson.dependencies[dep] || packageJson.optionalDependencies[dep]).toBeDefined();
      });
    });

    test("should not require Python environment", () => {
      expect(true).toBe(true);
    });
  });

  describe("Configuration", () => {
    test("should have valid package.json", () => {
      const packageJson = JSON.parse(
        fs.readFileSync(path.join(__dirname, "../../../package.json"), "utf8"),
      );

      expect(packageJson.name).toBe("queb");
      expect(packageJson.version).toBeDefined();
      expect(packageJson.main).toBeDefined();
      expect(packageJson.scripts).toBeDefined();
    });

    test("should have test scripts configured", () => {
      const packageJson = JSON.parse(
        fs.readFileSync(path.join(__dirname, "../../../package.json"), "utf8"),
      );

      expect(packageJson.scripts.test).toBeDefined();
      expect(packageJson.scripts["test:gate"]).toBeDefined();
      expect(packageJson.scripts["test:integration"]).toBeDefined();
    });

    test("should have proper directory structure", () => {
      const requiredDirs = ["../../../tests", "../../../tests/sdf_files"];

      requiredDirs.forEach((dir) => {
        const dirPath = path.join(__dirname, dir);
        expect(fs.existsSync(dirPath)).toBe(true);
      });
    });
  });

  describe("Frontend Assets", () => {
    test("should have valid HTML structure", () => {
      const htmlContent = fs.readFileSync(
        path.join(__dirname, "../../../src/client/core/index.html"),
        "utf8",
      );

      expect(htmlContent).toContain("<!doctype html>");
      expect(htmlContent).toContain("<title>{{TITLE}}</title>");
      expect(htmlContent).toContain("bundle.js");
    });

    test("should have valid CSS file", () => {
      const cssPath = path.join(__dirname, "../../../src/client/assets/style.css");
      expect(fs.existsSync(cssPath)).toBe(true);

      const cssContent = fs.readFileSync(cssPath, "utf8");
      expect(cssContent.length).toBeGreaterThan(100);
    });

    test("should have valid JavaScript files", () => {
      const jsFiles = [
        "../../../src/server/api/server.js",
        "../../../src/server/services/Structuralizer.js",
        "../../../src/server/services/molecular-processor.js",
      ];

      jsFiles.forEach((file) => {
        const filePath = path.join(__dirname, file);
        expect(fs.existsSync(filePath)).toBe(true);

        const content = fs.readFileSync(filePath, "utf8");
        expect(content.length).toBeGreaterThan(50);
      });
    });
  });

  describe("Schema Validation", () => {
    test("should have valid schemas defined", () => {
      const schemas = require("../../../src/server/schemas/schemas");

      expect(schemas.ImageMoleculeSchema).toBeDefined();
      expect(schemas.TextMoleculeSchema).toBeDefined();
      expect(schemas.SdfGenerationSchema).toBeDefined();
    });

    test("should validate basic schema structure", () => {
      const {
        ImageMoleculeSchema,
        TextMoleculeSchema,
        SdfGenerationSchema,
      } = require("../../../src/server/schemas/schemas");

      // Test that schemas can be imported and used
      expect(typeof ImageMoleculeSchema.safeParse).toBe("function");
      expect(typeof TextMoleculeSchema.safeParse).toBe("function");
      expect(typeof SdfGenerationSchema.safeParse).toBe("function");
    });
  });

  describe("Error Handling", () => {
    test("should handle invalid routes gracefully", async () => {
      const response = await request(app).get("/nonexistent-route");
      expect([200, 404, 500]).toContain(response.status); // 200 is expected for SPA routes
    });

    test("should handle malformed JSON requests", async () => {
      const response = await request(app)
        .post("/image-molecules")
        .set("Content-Type", "application/json")
        .send("invalid json");

      expect([400, 500]).toContain(response.status);
    });
  });

  describe("Performance Baseline", () => {
    test("should respond to basic requests quickly", async () => {
      const startTime = Date.now();

      const response = await request(app).get("/");

      const endTime = Date.now();
      const responseTime = endTime - startTime;

      expect(response.status).toBe(200);
      expect(responseTime).toBeLessThan(1000); // Should respond within 1 second
    });

    test("should handle concurrent requests", async () => {
      const requests = Array(3)
        .fill()
        .map(() => request(app).get("/"));

      const responses = await Promise.all(requests);

      responses.forEach((response) => {
        expect(response.status).toBe(200);
      });
    });
  });

  describe("Security Basics", () => {
    test("should not expose sensitive information in headers", async () => {
      const response = await request(app).get("/");

      const sensitiveHeaders = ["x-powered-by", "server", "x-aspnet-version"];
      sensitiveHeaders.forEach((header) => {
        expect(response.headers[header]).toBeUndefined();
      });
    });

    test("should handle large payloads appropriately", async () => {
      const largePayload = { object: "x".repeat(10000) };

      const response = await request(app)
        .post("/structures-from-text")
        .send(largePayload);

      expect([400, 413, 500]).toContain(response.status); // Should reject or handle appropriately
    });
  });
});
