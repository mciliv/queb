// test/fixtures/utils.js - Test utilities and helper functions
const fs = require("fs");
const path = require("path");
const { TEST_CONFIG, TEST_ENDPOINTS, TEST_FIXTURES } = require("./config");

// ==================== TEST UTILITIES ====================

/**
 * Sleep utility for tests
 * @param {number} ms - milliseconds to sleep
 */
function sleep(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

/**
 * Generate random test data
 */
function generateTestData() {
  return {
    randomSmiles: () => {
      const smiles = TEST_FIXTURES.molecules.simple;
      return smiles[Math.floor(Math.random() * smiles.length)];
    },

    randomObject: () => {
      const objects = TEST_FIXTURES.objects.common;
      return objects[Math.floor(Math.random() * objects.length)];
    },

    randomId: () => {
      return Math.random().toString(36).substr(2, 9);
    },

    timestamp: () => {
      return new Date().toISOString();
    },
  };
}

/**
 * Mock HTTP client for testing
 */
class MockHttpClient {
  constructor(baseUrl = TEST_ENDPOINTS.base) {
    this.baseUrl = baseUrl;
    this.responses = new Map();
    this.requests = [];
  }

  // Mock a response for a specific endpoint
  mockResponse(endpoint, response, status = 200) {
    this.responses.set(endpoint, { response, status });
  }

  // Clear all mocked responses
  clearMocks() {
    this.responses.clear();
    this.requests = [];
  }

  // Get request history
  getRequests() {
    return [...this.requests];
  }

  // Mock fetch implementation
  async fetch(endpoint, options = {}) {
    const url = `${this.baseUrl}${endpoint}`;
    const request = {
      url,
      endpoint,
      method: options.method || "GET",
      body: options.body,
      headers: options.headers,
      timestamp: new Date().toISOString(),
    };

    this.requests.push(request);

    // Check if we have a mocked response
    if (this.responses.has(endpoint)) {
      const mock = this.responses.get(endpoint);
      await sleep(TEST_CONFIG.mock.openai.delay);

      return {
        ok: mock.status >= 200 && mock.status < 300,
        status: mock.status,
        json: async () => mock.response,
        text: async () => JSON.stringify(mock.response),
      };
    }

    // If no mock, throw error
    throw new Error(`No mock response configured for ${endpoint}`);
  }
}

/**
 * Test file manager
 */
class TestFileManager {
  constructor() {
    this.tempFiles = [];
    this.testDir = TEST_CONFIG.directories.temp;
  }

  // Create a temporary test file
  createTempFile(filename, content) {
    // Ensure test directory exists
    if (!fs.existsSync(this.testDir)) {
      fs.mkdirSync(this.testDir, { recursive: true });
    }

    const filepath = path.join(this.testDir, filename);
    fs.writeFileSync(filepath, content);
    this.tempFiles.push(filepath);
    return filepath;
  }

  // Create a test SDF file
  createTestSdf(smiles, filename) {
    const sdfContent = `
  Mrv2114 12152024

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <SMILES>
${smiles}

> <_test_file>
true

$$$$
`;
    return this.createTempFile(
      filename || `${smiles.replace(/[^a-zA-Z0-9]/g, "_")}.sdf`,
      sdfContent,
    );
  }

  // Create a test image file
  createTestImage(imageData, filename) {
    const buffer = Buffer.from(imageData, "base64");
    const filepath = path.join(this.testDir, filename);
    fs.writeFileSync(filepath, buffer);
    this.tempFiles.push(filepath);
    return filepath;
  }

  // Clean up all temporary files
  cleanup() {
    this.tempFiles.forEach((filepath) => {
      if (fs.existsSync(filepath)) {
        fs.unlinkSync(filepath);
      }
    });
    this.tempFiles = [];
  }
}

/**
 * Test assertion utilities
 */
class TestAssertions {
  static isValidSmiles(smiles) {
    // Basic SMILES validation
    if (typeof smiles !== "string" || smiles.length === 0) {
      return false;
    }

    // SMILES should not contain spaces or certain invalid strings
    if (smiles.includes(" ") || smiles === "invalid") {
      return false;
    }

    // Basic SMILES character validation
    return /^[A-Za-z0-9\[\]()=#+\-\.@:\/\\%]+$/.test(smiles);
  }

  static isValidSdfPath(path) {
    return (
      typeof path === "string" &&
      path.endsWith(".sdf") &&
      (path.startsWith("/") || path.startsWith("http"))
    );
  }

  static hasValidTestResponse(response) {
    return (
      response &&
      response.output &&
      response.output._test &&
      response.output._test.timestamp
    );
  }

  static arrayContainsValidSmiles(array) {
    return (
      Array.isArray(array) &&
      array.length > 0 &&
      array.every((smiles) => this.isValidSmiles(smiles))
    );
  }

  static responseMatchesFixture(response, fixture) {
    if (!response || !response.output) return false;

    const { smiles } = response.output;
    const expectedSmiles = fixture.smiles;

    // Check if response contains expected SMILES
    return expectedSmiles.some((expected) => smiles.includes(expected));
  }
}

/**
 * Test data builder
 */
class TestDataBuilder {
  constructor() {
    this.data = {};
  }

  withObject(object) {
    this.data.object = object;
    return this;
  }

  withSmiles(smiles) {
    this.data.smiles = Array.isArray(smiles) ? smiles : [smiles];
    return this;
  }

  withImage(imageBase64) {
    this.data.imageBase64 = imageBase64;
    return this;
  }

  withCroppedImage(croppedImageBase64) {
    this.data.croppedImageBase64 = croppedImageBase64;
    return this;
  }

  withCoordinates(x, y) {
    this.data.x = x;
    this.data.y = y;
    return this;
  }

  withTestObject(testObject) {
    this.data.testObject = testObject;
    return this;
  }

  build() {
    return { ...this.data };
  }
}

/**
 * Test server utilities
 */
class TestServerUtils {
  static async waitForServer(port, timeout = 5000) {
    const start = Date.now();

    while (Date.now() - start < timeout) {
      try {
        const response = await fetch(`http://localhost:${port}/test/health`);
        if (response.ok) {
          return true;
        }
      } catch (error) {
        // Server not ready yet
      }

      await sleep(100);
    }

    return false;
  }

  static async resetTestEnvironment(port = TEST_CONFIG.server.port) {
    try {
      const response = await fetch(`http://localhost:${port}/test/utils/reset`);
      return response.ok;
    } catch (error) {
      console.error("Failed to reset test environment:", error);
      return false;
    }
  }

  static async getTestFixtures(port = TEST_CONFIG.server.port) {
    try {
      const response = await fetch(`http://localhost:${port}/test/fixtures`);
      return response.ok ? await response.json() : null;
    } catch (error) {
      console.error("Failed to get test fixtures:", error);
      return null;
    }
  }
}

/**
 * Test reporter
 */
class TestReporter {
  constructor() {
    this.results = [];
    this.startTime = Date.now();
  }

  logTest(testName, status, duration, details = {}) {
    const result = {
      testName,
      status,
      duration,
      details,
      timestamp: new Date().toISOString(),
    };

    this.results.push(result);
    console.log(`[TEST] ${testName}: ${status} (${duration}ms)`);
  }

  generateReport() {
    const endTime = Date.now();
    const totalDuration = endTime - this.startTime;

    const passed = this.results.filter((r) => r.status === "PASSED").length;
    const failed = this.results.filter((r) => r.status === "FAILED").length;
    const total = this.results.length;

    return {
      summary: {
        total,
        passed,
        failed,
        duration: totalDuration,
        passRate: total > 0 ? ((passed / total) * 100).toFixed(1) : "0",
      },
      results: this.results,
    };
  }
}

// ==================== EXPORTS ====================
module.exports = {
  sleep,
  generateTestData,
  MockHttpClient,
  TestFileManager,
  TestAssertions,
  TestDataBuilder,
  TestServerUtils,
  TestReporter,
};
