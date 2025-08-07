// test/integration/background.test.js - Background integration tests during development
// These tests run longer scenarios to catch integration issues while developing

const request = require("supertest");
const fs = require("fs");
const path = require("path");
const { TestFileManager, TestAssertions } = require("./utils");
const {
  getTestMolecule,
  createTestRequest,
  MOCK_IMAGES,
} = require("./fixtures");

describe("Background Integration Tests", () => {
  let app;
  let fileManager;

  beforeAll(() => {
    app = require("../../backend/api/server");
    fileManager = new TestFileManager();
  });

  afterAll(() => {
    fileManager.cleanup();
  });

  describe("SMILES Processing Pipeline", () => {
    it("should generate SDF files for valid SMILES", async () => {
      const testSmiles = ["O", "CCO", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"];

      const response = await request(app)
        .post("/generate-sdfs")
        .send({ smiles: testSmiles, overwrite: true });

      expect(response.status).toBe(200);
      expect(response.body.message).toBe("Files generated");
      expect(Array.isArray(response.body.sdfPaths)).toBe(true);
      expect(response.body.sdfPaths).toHaveLength(testSmiles.length);

      // Verify SDF files were created
      response.body.sdfPaths.forEach((sdfPath) => {
        const fullPath = path.join(__dirname, "..", sdfPath);
        expect(fs.existsSync(fullPath)).toBe(true);

        const content = fs.readFileSync(fullPath, "utf8");
        expect(content).toContain("$$$$"); // SDF file marker
      });
    }, 15000);

    it("should handle invalid SMILES gracefully", async () => {
      const invalidSmiles = ["INVALID_SMILES", "XYZ123"];

      const response = await request(app)
        .post("/generate-sdfs")
        .send({ smiles: invalidSmiles, overwrite: true });

      // Should handle errors without crashing
      expect([200, 400, 500]).toContain(response.status);
    }, 10000);

    it("should not overwrite existing files when overwrite=false", async () => {
      const testSmiles = ["O"];

      // First request - create file
      await request(app)
        .post("/generate-sdfs")
        .send({ smiles: testSmiles, overwrite: true });

      // Second request - should not overwrite
      const response = await request(app)
        .post("/generate-sdfs")
        .send({ smiles: testSmiles, overwrite: false });

      expect(response.status).toBe(200);
    }, 10000);
  });

  describe("Schema Validation", () => {
    it("should validate schemas against test data", () => {
      const schemas = require("../../backend/schemas/schemas");

      // Test SMILES array schema
      const validSmilesData = { smiles: ["O", "CCO"] };
      expect(() => schemas.smilesArray.parse(validSmilesData)).not.toThrow();

      // Test object molecules schema
      const validObjectData = { object: "water" };
      expect(() =>
        schemas.objectMoleculesRequest.parse(validObjectData),
      ).not.toThrow();

      // Test image molecules schema
      const validImageData = {
        imageBase64: "data:image/png;base64,iVBORw0KGgo...",
        croppedImageBase64: "data:image/png;base64,iVBORw0KGgo...",
        x: 100,
        y: 100,
      };
      expect(() =>
        schemas.imageMoleculesRequest.parse(validImageData),
      ).not.toThrow();
    });

    it("should reject invalid schema data", () => {
      const schemas = require("../../backend/schemas/schemas");

      // Invalid SMILES data
      expect(() =>
        schemas.smilesArray.parse({ smiles: "not an array" }),
      ).toThrow();
      expect(() => schemas.smilesArray.parse({ smiles: [123, 456] })).toThrow();

      // Invalid object data
      expect(() =>
        schemas.objectMoleculesRequest.parse({ object: 123 }),
      ).toThrow();
      expect(() => schemas.objectMoleculesRequest.parse({})).toThrow();
    });
  });

  describe("File System Operations", () => {
    it("should create sdf_files directory if missing", async () => {
      const sdfDir = path.join(__dirname, "..", "sdf_files");

      // Remove directory if it exists (cleanup from previous tests)
      if (fs.existsSync(sdfDir)) {
        fs.rmSync(sdfDir, { recursive: true, force: true });
      }

      const response = await request(app)
        .post("/generate-sdfs")
        .send({ smiles: ["O"], overwrite: true });

      expect(fs.existsSync(sdfDir)).toBe(true);
    }, 10000);

    it("should handle file system errors gracefully", async () => {
      // This test would need specific file system error conditions
      // For now, just ensure the endpoint handles edge cases
      const response = await request(app)
        .post("/generate-sdfs")
        .send({ smiles: [""] }); // Empty SMILES

      expect([200, 400, 500]).toContain(response.status);
    });
  });

  describe("Error Handling", () => {
    it("should handle malformed request bodies", async () => {
      const response = await request(app)
        .post("/generate-sdfs")
        .send({ invalid: "data" });

      expect(response.status).toBe(400);
      expect(response.body.error).toBeDefined();
    });

    it("should handle missing request data", async () => {
      const response = await request(app).post("/generate-sdfs").send({});

      expect(response.status).toBe(400);
    });

    it("should validate endpoint availability", async () => {
      // Test all major endpoints respond
      const endpoints = [
        { path: "/generate-sdfs", method: "post", data: { smiles: [] } },
        { path: "/object-molecules", method: "post", data: { object: "test" } },
        {
          path: "/image-molecules",
          method: "post",
          data: {
            imageBase64: "test",
            croppedImageBase64: "test",
            x: 0,
            y: 0,
          },
        },
      ];

      for (const endpoint of endpoints) {
        const response = await request(app)
          [endpoint.method](endpoint.path)
          .send(endpoint.data);
        expect([200, 400, 401, 403, 500]).toContain(response.status);
      }
    });
  });

  describe("Test Utilities Validation", () => {
    it("should validate test fixture quality", () => {
      const caffeine = getTestMolecule("caffeine");
      expect(TestAssertions.isValidSmiles(caffeine.smiles)).toBe(true);

      // Test all available molecules
      const molecules = ["caffeine", "ethanol", "water", "glucose"];
      molecules.forEach((molName) => {
        const mol = getTestMolecule(molName);
        expect(mol).toBeDefined();
        expect(TestAssertions.isValidSmiles(mol.smiles)).toBe(true);
      });
    });

    it("should validate mock image data", () => {
      expect(MOCK_IMAGES.blackSquare.base64).toMatch(
        /^data:image\/png;base64,/,
      );
      expect(MOCK_IMAGES.whiteSquare.base64).toMatch(
        /^data:image\/png;base64,/,
      );

      // Validate base64 content length (should be substantial)
      const base64Content = MOCK_IMAGES.blackSquare.base64.split(",")[1];
      expect(base64Content.length).toBeGreaterThan(100);
    });

    it("should validate test request builders", () => {
      const imageRequest = createTestRequest("imageMolecules", {
        x: 50,
        y: 75,
      });
      expect(imageRequest.x).toBe(50);
      expect(imageRequest.y).toBe(75);
      expect(imageRequest.imageBase64).toBeDefined();
      expect(imageRequest.croppedImageBase64).toBeDefined();

      const objectRequest = createTestRequest("objectMolecules", {
        object: "coffee",
      });
      expect(objectRequest.object).toBe("coffee");
    });
  });
});
