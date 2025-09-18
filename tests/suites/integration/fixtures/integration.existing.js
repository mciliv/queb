// test/fixtures/integration.existing.js - Integration tests using separated test functionality

const request = require("supertest");
const app = require("../../backend/api/server");
const {
  TEST_MOLECULES,
  TEST_OBJECTS,
  MOCK_IMAGES,
  createTestRequest,
  validateTestResponse,
  getTestMolecule,
  getTestObject,
} = require("./fixtures");
const { TestFileManager, TestAssertions, TestDataBuilder } = require("./utils");

describe("Separated Test Functionality", () => {
  let fileManager;

  beforeEach(() => {
    fileManager = new TestFileManager();
  });

  afterEach(() => {
    fileManager.cleanup();
    jest.clearAllMocks();
  });

  describe("Test Fixtures", () => {
    it("should provide valid test molecules", () => {
      const caffeine = getTestMolecule("caffeine");

      expect(caffeine).toBeDefined();
      expect(caffeine.smiles).toBe("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
      expect(caffeine.name).toBe("Caffeine");
      expect(TestAssertions.isValidSmiles(caffeine.smiles)).toBe(true);
    });

    it("should provide valid test objects", () => {
      const coffee = getTestObject("coffee");

      expect(coffee).toBeDefined();
      expect(coffee.object).toBe("coffee");
      expect(Array.isArray(coffee.chemicals)).toBe(true);

      // Extract SMILES from chemicals for testing
      const smiles = coffee.chemicals
        .map((chem) => chem.smiles)
        .filter(Boolean);
      expect(Array.isArray(smiles)).toBe(true);
      expect(TestAssertions.arrayContainsValidSmiles(smiles)).toBe(true);
    });

    it("should provide mock images", () => {
      expect(MOCK_IMAGES.blackSquare.base64).toBeDefined();
      expect(MOCK_IMAGES.whiteSquare.base64).toBeDefined();
      expect(typeof MOCK_IMAGES.blackSquare.base64).toBe("string");
    });
  });

  describe("Test Data Builder", () => {
    it("should build test request data", () => {
      const testData = new TestDataBuilder()
        .withObject("coffee")
        .withSmiles(["CN1C=NC2=C1C(=O)N(C(=O)N2C)C"])
        .withCoordinates(100, 200)
        .build();

      expect(testData.object).toBe("coffee");
      expect(testData.smiles).toEqual(["CN1C=NC2=C1C(=O)N(C(=O)N2C)C"]);
      expect(testData.x).toBe(100);
      expect(testData.y).toBe(200);
    });
  });

  describe("Test Request Creation", () => {
    it("should create image molecules test request", () => {
      const request = createTestRequest("imageMolecules", {
        object: "coffee",
        x: 150,
        y: 250,
      });

      expect(request.imageBase64).toBeDefined();
      expect(request.croppedImageBase64).toBeDefined();
      expect(request.x).toBe(150);
      expect(request.y).toBe(250);
    });

    it("should create object molecules test request", () => {
      const request = createTestRequest("objectMolecules", {
        object: "wine",
      });

      expect(request.object).toBe("wine");
    });
  });

  describe("Response Validation", () => {
    it("should validate image molecules response", () => {
      const validResponse = {
        output: {
          object: "coffee",
          chemicals: [
            { name: "Caffeine", smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" },
            { name: "Water", smiles: "O" },
          ],
        },
      };

      expect(validateTestResponse(validResponse, "imageMolecules")).toBe(true);
    });

    it("should validate object molecules response", () => {
      const validResponse = {
        output: {
          object: "wine",
          chemicals: [
            { name: "Ethanol", smiles: "CCO" },
            { name: "Water", smiles: "O" },
          ],
        },
      };

      expect(validateTestResponse(validResponse, "objectMolecules")).toBe(true);
    });

    it("should reject invalid responses", () => {
      const invalidResponse = {
        output: {
          object: "coffee",
          // missing chemicals array
        },
      };

      expect(validateTestResponse(invalidResponse, "imageMolecules")).toBe(
        false,
      );
    });
  });

  describe("File Manager", () => {
    it("should create and cleanup test files", () => {
      const fs = require("fs");
      const testContent = "test content";

      const filepath = fileManager.createTempFile("test.txt", testContent);

      expect(fs.existsSync(filepath)).toBe(true);
      expect(fs.readFileSync(filepath, "utf8")).toBe(testContent);

      fileManager.cleanup();
      expect(fs.existsSync(filepath)).toBe(false);
    });

    it("should create test SDF files", () => {
      const fs = require("fs");
      const smiles = "CCO";

      const sdfPath = fileManager.createTestSdf(smiles, "ethanol.sdf");

      expect(fs.existsSync(sdfPath)).toBe(true);

      const content = fs.readFileSync(sdfPath, "utf8");
      expect(content).toContain(smiles);
      expect(content).toContain("_test_file");

      fileManager.cleanup();
      expect(fs.existsSync(sdfPath)).toBe(false);
    });
  });

  describe("Test Assertions", () => {
    it("should validate SMILES strings", () => {
      expect(TestAssertions.isValidSmiles("CCO")).toBe(true);
      expect(TestAssertions.isValidSmiles("O")).toBe(true);
      expect(TestAssertions.isValidSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")).toBe(
        true,
      );

      expect(TestAssertions.isValidSmiles("")).toBe(false);
      expect(TestAssertions.isValidSmiles("invalid smiles")).toBe(false);
      expect(TestAssertions.isValidSmiles(123)).toBe(false);
    });

    it("should validate SDF paths", () => {
      expect(TestAssertions.isValidSdfPath("/path/to/file.sdf")).toBe(true);
      expect(TestAssertions.isValidSdfPath("http://example.com/file.sdf")).toBe(
        true,
      );

      expect(TestAssertions.isValidSdfPath("file.txt")).toBe(false);
      expect(TestAssertions.isValidSdfPath("")).toBe(false);
      expect(TestAssertions.isValidSdfPath(null)).toBe(false);
    });

    it("should validate arrays of SMILES", () => {
      expect(TestAssertions.arrayContainsValidSmiles(["CCO", "O"])).toBe(true);

      expect(TestAssertions.arrayContainsValidSmiles([])).toBe(false);
      expect(TestAssertions.arrayContainsValidSmiles(["CCO", "invalid"])).toBe(
        false,
      );
      expect(TestAssertions.arrayContainsValidSmiles("not an array")).toBe(
        false,
      );
    });
  });

  describe("Integration with existing endpoints", () => {
    // Skip these tests for now due to OpenAI mocking complexity
    // These can be enabled once proper mocking is set up
    it.skip("should work with text analysis endpoint", async () => {
      const testRequest = createTestRequest("objectMolecules", {
        object: "coffee",
      });

      const response = await request(app)
        .post("/object-molecules")
        .send(testRequest);

      expect(response.status).toBe(200);
      expect(validateTestResponse(response.body, "objectMolecules")).toBe(true);
    });

    it.skip("should work with image analysis endpoint", async () => {
      const testRequest = createTestRequest("imageMolecules", {
        imageBase64: MOCK_IMAGES.blackSquare.base64,
        croppedImageBase64: MOCK_IMAGES.whiteSquare.base64,
        x: 100,
        y: 100,
      });

      const response = await request(app)
        .post("/image-molecules")
        .send(testRequest);

      expect(response.status).toBe(200);
      expect(validateTestResponse(response.body, "imageMolecules")).toBe(true);
    });

    it("should validate that test fixtures work for future endpoint integration", () => {
      // Test that our fixtures and utilities are ready for integration testing
      const textRequest = createTestRequest("objectMolecules", {
        object: "coffee",
      });
      const imageRequest = createTestRequest("imageMolecules");

      expect(textRequest.object).toBe("coffee");
      expect(imageRequest.imageBase64).toBeDefined();
      expect(imageRequest.x).toBeDefined();
      expect(imageRequest.y).toBeDefined();

      // Mock response validation
      const mockResponse = {
        output: {
          object: "test object",
          chemicals: [
            { name: "Water", smiles: "O" },
            { name: "Ethanol", smiles: "CCO" },
          ],
        },
      };

      expect(validateTestResponse(mockResponse, "imageMolecules")).toBe(true);
    });
  });
});
