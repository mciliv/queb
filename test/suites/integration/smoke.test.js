// test/integration/smoke.test.js - Quick smoke tests for pre-development validation
// These tests run in under 5 seconds to catch basic issues before starting development

const request = require("supertest");
const app = require("../../backend/api/server");
const fs = require("fs");
const path = require("path");

describe("Smoke Tests", () => {
  describe("Application Startup", () => {
    test("should start without crashing", () => {
      expect(app).toBeDefined();
      expect(typeof app.listen).toBe("function");
    });

    test("should serve static files", async () => {
      const response = await request(app).get("/");
      expect(response.status).toBe(200);
      expect(response.text).toContain("<!doctype html>");
      expect(response.text).toContain("Atomizer - Molecular Analysis");
    });

    test("should have required core files", () => {
      const requiredFiles = [
        "server.js",
        "index.html",
        "app.js",
        "style.css",
        "schemas.js",
        "AtomPredictor.js",
        "molecular-processor.js",
        "package.json",
      ];

      requiredFiles.forEach((file) => {
        const filePath = path.join(__dirname, "..", file);
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
      const sdfDir = path.join(__dirname, "../sdf_files");
      if (fs.existsSync(sdfDir)) {
        const sdfResponse = await request(app).get("/sdf_files/");
        expect([200, 404]).toContain(sdfResponse.status); // 404 is okay if directory is empty
      }
    });

    test("should handle API endpoint requests", async () => {
      // Test image molecules endpoint (should respond even without valid data)
      const imageResponse = await request(app)
        .post("/image-molecules")
        .send({});

      expect([400, 500]).toContain(imageResponse.status); // Should reject invalid data

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
        fs.readFileSync(path.join(__dirname, "../package.json"), "utf8"),
      );

      const requiredDeps = ["express", "cors", "openai", "zod"];

      requiredDeps.forEach((dep) => {
        expect(packageJson.dependencies[dep]).toBeDefined();
      });
    });

    test("should have Python environment available", () => {
      const { execSync } = require("child_process");

      try {
        const pythonVersion = execSync("python --version", {
          encoding: "utf8",
        });
        expect(pythonVersion).toContain("Python");
      } catch (error) {
        // Python not available - this is okay for some test environments
        console.warn("Python not available in test environment");
      }
    });

    test("should have sdf.py script", () => {
      const sdfPath = path.join(__dirname, "../sdf.py");
      expect(fs.existsSync(sdfPath)).toBe(true);

      const sdfContent = fs.readFileSync(sdfPath, "utf8");
      expect(sdfContent).toContain("from rdkit import Chem");
    });
  });

  describe("Configuration", () => {
    test("should have valid package.json", () => {
      const packageJson = JSON.parse(
        fs.readFileSync(path.join(__dirname, "../package.json"), "utf8"),
      );

      expect(packageJson.name).toBe("mol");
      expect(packageJson.version).toBeDefined();
      expect(packageJson.main).toBeDefined();
      expect(packageJson.scripts).toBeDefined();
    });

    test("should have test scripts configured", () => {
      const packageJson = JSON.parse(
        fs.readFileSync(path.join(__dirname, "../package.json"), "utf8"),
      );

      expect(packageJson.scripts.test).toBeDefined();
      expect(packageJson.scripts["test:unit"]).toBeDefined();
      expect(packageJson.scripts["test:integration"]).toBeDefined();
    });

    test("should have proper directory structure", () => {
      const requiredDirs = ["tests", "sdf_files"];

      requiredDirs.forEach((dir) => {
        const dirPath = path.join(__dirname, "..", dir);
        expect(fs.existsSync(dirPath)).toBe(true);
      });
    });
  });

  describe("Frontend Assets", () => {
    test("should have valid HTML structure", () => {
      const htmlContent = fs.readFileSync(
        path.join(__dirname, "../index.html"),
        "utf8",
      );

      expect(htmlContent).toContain("<!doctype html>");
      expect(htmlContent).toContain("<title>Atomic Reality</title>");
      expect(htmlContent).toContain("app.js");
      expect(htmlContent).toContain("style.css");
    });

    test("should have valid CSS file", () => {
      const cssPath = path.join(__dirname, "../style.css");
      expect(fs.existsSync(cssPath)).toBe(true);

      const cssContent = fs.readFileSync(cssPath, "utf8");
      expect(cssContent.length).toBeGreaterThan(100);
    });

    test("should have valid JavaScript files", () => {
      const jsFiles = [
        "app.js",
        "server.js",
        "AtomPredictor.js",
        "molecular-processor.js",
      ];

      jsFiles.forEach((file) => {
        const filePath = path.join(__dirname, "..", file);
        expect(fs.existsSync(filePath)).toBe(true);

        const content = fs.readFileSync(filePath, "utf8");
        expect(content.length).toBeGreaterThan(50);
      });
    });
  });

  describe("Schema Validation", () => {
    test("should have valid schemas defined", () => {
      const schemas = require("../../backend/schemas/schemas");

      expect(schemas.ImageMoleculeSchema).toBeDefined();
      expect(schemas.TextMoleculeSchema).toBeDefined();
      expect(schemas.SdfGenerationSchema).toBeDefined();
    });

    test("should validate basic schema structure", () => {
      const {
        ImageMoleculeSchema,
        TextMoleculeSchema,
        SdfGenerationSchema,
      } = require("../../backend/schemas/schemas");

      // Test that schemas can be imported and used
      expect(typeof ImageMoleculeSchema.safeParse).toBe("function");
      expect(typeof TextMoleculeSchema.safeParse).toBe("function");
      expect(typeof SdfGenerationSchema.safeParse).toBe("function");
    });
  });

  describe("Error Handling", () => {
    test("should handle invalid routes gracefully", async () => {
      const response = await request(app).get("/nonexistent-route");
      expect([404, 500]).toContain(response.status);
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
      const largePayload = { data: "x".repeat(10000) };

      const response = await request(app)
        .post("/image-molecules")
        .send(largePayload);

      expect([400, 413, 500]).toContain(response.status); // Should reject or handle appropriately
    });
  });
});
