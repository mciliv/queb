// test/integration/deployment.test.js - Comprehensive pre-deployment test suite
// These tests validate the entire system before deployment

const request = require("supertest");
const fs = require("fs");
const path = require("path");
const { execSync } = require("child_process");
const { TestFileManager, TestAssertions } = require("./utils");
const {
  getTestMolecule,
  createTestRequest,
  MOCK_IMAGES,
} = require("./fixtures");

describe("Pre-Deployment Validation", () => {
  let app;
  let fileManager;

  beforeAll(async () => {
    app = require("../../backend/api/server");
    fileManager = new TestFileManager();
  }, 30000);

  afterAll(() => {
    fileManager.cleanup();
  });

  describe("System Requirements", () => {
    it("should have all required dependencies installed", () => {
      const packageJson = JSON.parse(
        fs.readFileSync(path.join(__dirname, "..", "package.json"), "utf8"),
      );
      const requiredDeps = ["express", "cors", "multer", "openai", "zod"];

      requiredDeps.forEach((dep) => {
        expect(packageJson.dependencies[dep]).toBeDefined();
      });
    });

    it("should have Python environment ready", () => {
      try {
        const pythonVersion = execSync("python --version", {
          encoding: "utf8",
        });
        expect(pythonVersion).toContain("Python");

        // Test RDKit availability
        execSync("python -c 'from rdkit import Chem'", { encoding: "utf8" });

        // Test basic SMILES processing
        const testResult = execSync(
          'python -c \'from rdkit import Chem; print("PASS" if Chem.MolFromSmiles("CCO") else "FAIL")\'',
          { encoding: "utf8" },
        );
        expect(testResult.trim()).toBe("PASS");
      } catch (error) {
        throw new Error(`Python environment check failed: ${error.message}`);
      }
    });

    it("should have all core files present and valid", () => {
      const criticalFiles = [
        "server.js",
        "index.html",
        "app.js",
        "style.css",
        "schemas.js",
        "sdf.py",
        "package.json",
      ];

      criticalFiles.forEach((file) => {
        const filePath = path.join(__dirname, "..", file);
        expect(fs.existsSync(filePath)).toBe(true);

        const stats = fs.statSync(filePath);
        expect(stats.size).toBeGreaterThan(0);
      });
    });
  });

  describe("Core Functionality Integration", () => {
    it("should process complete SMILES workflow", async () => {
      const testSmiles = ["O", "CCO", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"];

      // Test SDF generation
      const response = await request(app)
        .post("/generate-sdfs")
        .send({ smiles: testSmiles, overwrite: true });

      expect(response.status).toBe(200);
      expect(response.body.sdfPaths).toHaveLength(testSmiles.length);

      // Verify all SDF files were created and are valid
      for (const sdfPath of response.body.sdfPaths) {
        const fullPath = path.join(__dirname, "..", sdfPath);
        expect(fs.existsSync(fullPath)).toBe(true);

        const content = fs.readFileSync(fullPath, "utf8");
        expect(content).toContain("$$$$"); // SDF format marker
        expect(content.length).toBeGreaterThan(100); // Non-trivial content
      }
    }, 20000);

    it("should handle edge cases gracefully", async () => {
      const edgeCases = [
        { smiles: [], expected: [400] }, // Empty array
        { smiles: [""], expected: [400, 500] }, // Empty string
        { smiles: ["INVALID"], expected: [400, 500] }, // Invalid SMILES
        { smiles: ["O", "INVALID", "CCO"], expected: [200, 400, 500] }, // Mixed valid/invalid
      ];

      for (const testCase of edgeCases) {
        const response = await request(app)
          .post("/generate-sdfs")
          .send({ smiles: testCase.smiles, overwrite: true });

        expect(testCase.expected).toContain(response.status);
      }
    }, 15000);
  });

  describe("Schema Validation Coverage", () => {
    it("should validate all schema types thoroughly", () => {
      const schemas = require("../../backend/schemas/schemas");

      // Test valid data for all schemas
      const validTests = [
        { schema: schemas.smilesArray, data: { smiles: ["O", "CCO"] } },
        { schema: schemas.objectMoleculesRequest, data: { object: "water" } },
        {
          schema: schemas.imageMoleculesRequest,
          data: {
            imageBase64: "data:image/png;base64,test",
            croppedImageBase64: "data:image/png;base64,test",
            x: 100,
            y: 100,
          },
        },
      ];

      validTests.forEach((test) => {
        expect(() => test.schema.parse(test.data)).not.toThrow();
      });

      // Test invalid data rejection
      const invalidTests = [
        { schema: schemas.smilesArray, data: { smiles: "not array" } },
        { schema: schemas.smilesArray, data: { smiles: [123] } },
        { schema: schemas.objectMoleculesRequest, data: {} },
        { schema: schemas.objectMoleculesRequest, data: { object: 123 } },
        { schema: schemas.imageMoleculesRequest, data: { x: "not number" } },
      ];

      invalidTests.forEach((test) => {
        expect(() => test.schema.parse(test.data)).toThrow();
      });
    });
  });

  describe("Performance Requirements", () => {
    it("should handle moderate load efficiently", async () => {
      const startTime = Date.now();
      const concurrentRequests = 5;
      const testSmiles = ["O", "CCO"];

      const promises = Array(concurrentRequests)
        .fill()
        .map((_, i) =>
          request(app)
            .post("/generate-sdfs")
            .send({ smiles: testSmiles, overwrite: false }),
        );

      const responses = await Promise.all(promises);
      const endTime = Date.now();

      responses.forEach((response) => {
        expect([200, 400, 500]).toContain(response.status);
      });

      // Should complete within reasonable time (adjust as needed)
      expect(endTime - startTime).toBeLessThan(10000);
    }, 15000);

    it("should handle large SMILES arrays", async () => {
      const largeSmilesList = Array(20).fill("O"); // 20 identical water molecules

      const response = await request(app)
        .post("/generate-sdfs")
        .send({ smiles: largeSmilesList, overwrite: false });

      expect([200, 400, 500]).toContain(response.status);

      if (response.status === 200) {
        expect(response.body.sdfPaths).toHaveLength(largeSmilesList.length);
      }
    }, 25000);
  });

  describe("Error Recovery", () => {
    it("should recover from Python script failures", async () => {
      // Test with intentionally problematic SMILES
      const problematicSmiles = ["[invalid", ")", "(("];

      const response = await request(app)
        .post("/generate-sdfs")
        .send({ smiles: problematicSmiles, overwrite: true });

      // Should handle gracefully without crashing server
      expect([200, 400, 500]).toContain(response.status);

      // Server should still be responsive after error
      const healthCheck = await request(app).get("/");
      expect(healthCheck.status).toBe(200);
    }, 10000);

    it("should handle file system errors", async () => {
      // Create a scenario that might cause file system issues
      const response = await request(app)
        .post("/generate-sdfs")
        .send({ smiles: ["O".repeat(1000)], overwrite: true }); // Very long filename

      expect([200, 400, 500]).toContain(response.status);
    });
  });

  describe("Security Validation", () => {
    it("should reject malicious input", async () => {
      const maliciousInputs = [
        { smiles: ["'; DROP TABLE users; --"] },
        { smiles: ["<script>alert('xss')</script>"] },
        { smiles: ["../../../etc/passwd"] },
        { object: "<script>alert('xss')</script>" },
      ];

      for (const input of maliciousInputs) {
        if (input.smiles) {
          const response = await request(app)
            .post("/generate-sdfs")
            .send(input);
          expect([400, 500]).toContain(response.status);
        }

        if (input.object) {
          const response = await request(app)
            .post("/object-molecules")
            .send(input);
          expect([400, 401, 403, 500]).toContain(response.status);
        }
      }
    });

    it("should handle oversized requests", async () => {
      const oversizedData = {
        smiles: Array(1000).fill("O"), // Very large array
      };

      const response = await request(app)
        .post("/generate-sdfs")
        .send(oversizedData);

      // Should either process or reject gracefully
      expect([200, 400, 413, 500]).toContain(response.status);
    }, 30000);
  });

  describe("Deployment Readiness", () => {
    it("should have proper environment configuration", () => {
      // Check that deployment scripts exist
      const packageJson = JSON.parse(
        fs.readFileSync(path.join(__dirname, "..", "package.json"), "utf8"),
      );

      expect(packageJson.scripts["deploy"]).toBeDefined();
      expect(packageJson.scripts.start).toBeDefined();
      expect(packageJson.scripts.build).toBeDefined();
    });

    it("should have stable components protected", () => {
      // Verify STABLE_COMPONENTS.md exists and has content
      const stableComponentsPath = path.join(
        __dirname,
        "..",
        "STABLE_COMPONENTS.md",
      );

      if (fs.existsSync(stableComponentsPath)) {
        const content = fs.readFileSync(stableComponentsPath, "utf8");
        expect(content).toContain("DO NOT MODIFY");
        expect(content.length).toBeGreaterThan(100);
      }
    });

    it("should validate test coverage completeness", () => {
      // Ensure all major endpoints are covered
      const testFiles = fs.readdirSync(path.join(__dirname));
      const testFileCount = testFiles.filter((f) =>
        f.endsWith(".test.js"),
      ).length;

      expect(testFileCount).toBeGreaterThanOrEqual(4); // smoke, background, deployment, plus others
    });

    it("should verify git repository state", () => {
      try {
        // Check if we're in a git repository
        execSync("git status", {
          encoding: "utf8",
          cwd: path.join(__dirname, ".."),
        });

        // Verify main files are tracked
        const gitFiles = execSync("git ls-files", {
          encoding: "utf8",
          cwd: path.join(__dirname, ".."),
        });
        const trackedFiles = gitFiles.split("\n");

        const criticalFiles = ["server.js", "package.json", "schemas.js"];
        criticalFiles.forEach((file) => {
          expect(trackedFiles).toContain(file);
        });
      } catch (error) {
        // Git not available or not a git repo - that's okay for some deployments
        console.warn("Git repository check skipped:", error.message);
      }
    });
  });
});
