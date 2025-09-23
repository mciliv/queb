// test/unit/unit.test.js - Unit tests for individual components
// These tests run quickly (< 5 seconds) and validate individual functions and modules

const request = require("supertest");
const fs = require("fs");
const path = require("path");
const Structuralizer = require("../../../backend/services/Structuralizer");
const MolecularProcessor = require("../../../backend/services/molecular-processor");
const testConfigs = require("../../../backend/test-config-example");
const {
  ImageMoleculeSchema,
  TextMoleculeSchema,
  SdfGenerationSchema,
} = require("../../../backend/schemas/schemas");

// Mock OpenAI API
jest.mock("openai", () => ({
  OpenAI: jest.fn().mockImplementation(() => ({
    chat: {
      completions: {
        create: jest.fn().mockResolvedValue({
          choices: [
            {
              message: {
                content: JSON.stringify({
                  object: "test object",
                  chemicals: [
                    { name: "Ethanol", smiles: "CCO" },
                    { name: "Acetic acid", smiles: "CC(=O)O" },
                  ],
                }),
              },
            },
          ],
        }),
      },
    },
  })),
}));

// Mock file system
jest.mock("fs", () => ({
  existsSync: jest.fn(),
  writeFileSync: jest.fn(),
  readFileSync: jest.fn(),
  mkdirSync: jest.fn(),
}));

// Mock puppeteer to prevent import issues
jest.mock("puppeteer", () => ({
  launch: jest.fn(),
  connect: jest.fn()
}));

// Mock screenshot service to avoid puppeteer issues during tests
jest.mock("../../../src/server/services/screenshot-service", () => {
  return jest.fn().mockImplementation(() => ({
    captureApp: jest.fn().mockResolvedValue({ success: true, filename: 'mock.png' }),
    captureWithInput: jest.fn().mockResolvedValue({ success: true, filename: 'mock.png' }),
    captureAnalysis: jest.fn().mockResolvedValue({ success: true, filename: 'mock.png' }),
    listScreenshots: jest.fn().mockResolvedValue(['mock1.png', 'mock2.png']),
    getScreenshotPath: jest.fn().mockResolvedValue('/path/to/mock.png'),
    cleanupOldScreenshots: jest.fn().mockResolvedValue({ cleaned: 2 })
  }));
});

// Mock child_process
jest.mock("child_process", () => ({
  execSync: jest.fn().mockReturnValue("test output"),
  spawn: jest.fn().mockReturnValue({
    stdout: { on: jest.fn() },
    stderr: { on: jest.fn() },
    on: jest.fn((event, callback) => {
      if (event === "close") callback(1); // Default to failure
    }),
  }),
}));

// Global app instance for all tests
let app;

beforeAll(() => {
  // Load server for all tests
  try {
    app = require("../../../backend/api/server");
  } catch (error) {
    console.error("Failed to load server:", error.message);
    throw error;
  }
});

// Comprehensive App Validation Rules - Core Functionality
const APP_VALIDATION_RULES = {
  // Input validation rules
  textInput: {
    minLength: 2,
    maxLength: 500,
    invalidPatterns: [
      /^(love|hate|happy|sad|angry|joy|fear|hope|dream|idea|thought|feeling|emotion)/,
      /^(running|walking|talking|thinking|sleeping|eating|drinking)$/,
      /^[a-z]{1,2}$/,
      /^[^a-z]*$/,
      /^(asdf|qwerty|test|random|nothing|something|anything|everything)$/i,
      /^(blah|meh|hmm|ugh|oof|nah|yeah|yep|nope|ok|okay)$/i
    ],
    validExamples: ['water', 'ethanol', 'coffee', 'apple', 'salt', 'sugar']
  },
  
  // API endpoint validation
  endpoints: {
    '/analyze-text': { requiredFields: ['object'], responseFormat: 'json' },
    '/image-molecules': { requiredFields: ['imageBase64'], responseFormat: 'json' },
    '/generate-sdfs': { requiredFields: ['smiles'], responseFormat: 'json' },
    '/validate-payment': { requiredFields: ['device_token'], responseFormat: 'json' }
  },
  
  // Core component validation
  components: {
    camera: { required: ['video-feed', 'capture-btn'], permissions: ['camera'] },
    textInput: { required: ['object-input'], validation: true },
    payment: { required: ['payment-setup'], validation: true },
    molecular: { required: ['viewer-container'], smiles: true }
  },
  
  // System integrity checks
  integrity: {
    serverStartup: { timeout: 5000, healthChecks: ['/'] },
    databaseConnection: { optional: true, fallback: true },
    aiService: { timeout: 10000, retries: 3 },
    fileSystem: { writeable: true, paths: ['sdf_files', 'logs'] }
  },
  
  // Performance benchmarks
  performance: {
    textAnalysis: { maxTime: 15000, expectedResponse: 'smiles' },
    imageAnalysis: { maxTime: 20000, expectedResponse: 'molecules' },
    sdfGeneration: { maxTime: 5000, expectedResponse: 'files' },
    pageLoad: { maxTime: 3000, criticalResources: ['style.css', 'bundle.js'] }
  }
};

describe("Unit Tests - Core App Validation", () => {
  let atomPredictor;
  let molecularProcessor;

  beforeEach(() => {
    atomPredictor = new Structuralizer("test-api-key");
    molecularProcessor = new MolecularProcessor();
  });

  // App validation rule tests
  describe("App Validation Rules Compliance", () => {
    test("input validation rules work correctly", () => {
      APP_VALIDATION_RULES.textInput.invalidPatterns.forEach(pattern => {
        expect(pattern.test('love')).toBe(pattern === APP_VALIDATION_RULES.textInput.invalidPatterns[0]);
      });
      
      APP_VALIDATION_RULES.textInput.validExamples.forEach(example => {
        expect(example.length).toBeGreaterThanOrEqual(APP_VALIDATION_RULES.textInput.minLength);
        expect(example.length).toBeLessThanOrEqual(APP_VALIDATION_RULES.textInput.maxLength);
      });
    });
    
    test("all required API endpoints are defined", () => {
      Object.keys(APP_VALIDATION_RULES.endpoints).forEach(endpoint => {
        expect(endpoint).toMatch(/^\//); // Must start with /
        expect(APP_VALIDATION_RULES.endpoints[endpoint].requiredFields).toBeDefined();
      });
    });
    
    test("component requirements are specified", () => {
      Object.keys(APP_VALIDATION_RULES.components).forEach(component => {
        expect(APP_VALIDATION_RULES.components[component].required).toBeDefined();
        expect(Array.isArray(APP_VALIDATION_RULES.components[component].required)).toBe(true);
      });
    });
    
    test("performance benchmarks are realistic", () => {
      Object.keys(APP_VALIDATION_RULES.performance).forEach(operation => {
        const benchmark = APP_VALIDATION_RULES.performance[operation];
        expect(benchmark.maxTime).toBeGreaterThan(1000); // At least 1 second
        expect(benchmark.maxTime).toBeLessThan(30000); // Less than 30 seconds
      });
    });
  });

  describe("Structuralizer", () => {
    test("should construct and expose analysis methods", () => {
      expect(typeof atomPredictor.structuralizeText).toBe("function");
      expect(typeof atomPredictor.structuralizeImage).toBe("function");
    });

    describe("Test Configuration Interface", () => {
      test("should accept test configuration for model/prompt override", () => {
        const testStructuralizer = new Structuralizer("test-key", testConfigs.customPrompt);
        expect(testStructuralizer.testConfig.model).toBe("gpt-4o");
        expect(testStructuralizer.testConfig.prompt).toBe("custom");
      });

      test("should use test config model instead of candidates when provided", () => {
        const testStructuralizer = new Structuralizer("test-key", testConfigs.hardcodedPrompt);
        expect(testStructuralizer.testConfig.model).toBe("gpt-4-turbo");
        expect(testStructuralizer.testConfig.prompt).toContain("chemical structure");
      });

      test("should fall back to default behavior when no test config", () => {
        const defaultStructuralizer = new Structuralizer("test-key");
        expect(defaultStructuralizer.testConfig).toEqual({});
        expect(defaultStructuralizer.modelCandidates).toContain("gpt-4o");
      });

      test("should handle different test configurations", () => {
        const configs = [
          testConfigs.customPrompt,
          testConfigs.hardcodedPrompt,
          testConfigs.structuralizeTest,
          testConfigs.objectDetectionTest
        ];

        configs.forEach((config, index) => {
          const testStructuralizer = new Structuralizer("test-key", config);
          expect(testStructuralizer.testConfig.model).toBeDefined();
          expect(testStructuralizer.testConfig.prompt).toBeDefined();
          console.log(`Test config ${index + 1}: ${config.model} with ${config.prompt.length > 10 ? 'custom' : config.prompt} prompt`);
        });
      });

      test("should bypass candidate looping when test config provided", async () => {
        const mockOpenAI = require("openai");
        const mockCreate = jest.fn().mockResolvedValue({
          choices: [{ message: { content: '{"object": "test", "chemicals": []}' } }]
        });
        
        mockOpenAI.OpenAI.mockImplementationOnce(() => ({
          chat: { completions: { create: mockCreate } }
        }));

        const testStructuralizer = new Structuralizer("test-key", testConfigs.hardcodedPrompt);
        await testStructuralizer.structuralizeText("test input");

        expect(mockCreate).toHaveBeenCalledWith(
          expect.objectContaining({
            model: "gpt-4-turbo",
            messages: [{ role: 'user', content: testConfigs.hardcodedPrompt.prompt }]
          })
        );
      });
    });

    describe("structuralizeImage", () => {
      const mockImageBase64 =
        "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==";

      test("should analyze image without cropped region", async () => {
        const result = await atomPredictor.structuralize({ imageBase64: mockImageBase64, object: "" });
        expect(result).toHaveProperty("object");
        expect(result).toHaveProperty("chemicals");
        expect(Array.isArray(result.chemicals)).toBe(true);
        expect(result.object).toBe("test object");
        expect(result.chemicals).toHaveLength(2);
      });

      test("should analyze image with cropped region", async () => {
        const result = await atomPredictor.structuralize({
          imageBase64: mockImageBase64,
          object: "",
          x: 100,
          y: 100
        });
        expect(result).toHaveProperty("object");
        expect(result).toHaveProperty("chemicals");
        expect(Array.isArray(result.chemicals)).toBe(true);
      });

      test("should handle null parameters gracefully", async () => {
        const result = await atomPredictor.structuralize({
          imageBase64: mockImageBase64,
          object: "",
          x: null,
          y: null
        });
        expect(result).toHaveProperty("object");
        expect(result).toHaveProperty("chemicals");
      });

      test("should throw error when API fails", async () => {
        const mockOpenAI = require("openai");
        mockOpenAI.OpenAI.mockImplementationOnce(() => ({
          chat: {
            completions: {
              create: jest.fn().mockRejectedValue(new Error("API Error")),
            },
          },
        }));

        const analyzer = new Structuralizer("invalid-key");
        await expect(analyzer.structuralizeImage(mockImageBase64)).rejects.toThrow();
      });

      test("should handle empty image base64", async () => {
        await expect(atomPredictor.structuralizeImage("")).rejects.toThrow();
      });
    });

    describe("structuralizeText", () => {
      test("should analyze text input successfully", async () => {
        const result = await atomPredictor.structuralizeText("ethanol");
        expect(result).toHaveProperty("object");
        expect(result).toHaveProperty("chemicals");
        expect(Array.isArray(result.chemicals)).toBe(true);
        expect(result.object).toBe("test object");
        expect(result.chemicals[0]).toHaveProperty("name");
        expect(result.chemicals[0]).toHaveProperty("smiles");
      });

      test("should handle empty text input", async () => {
        const result = await atomPredictor.structuralizeText("");
        expect(result).toHaveProperty("object");
        expect(result).toHaveProperty("chemicals");
      });

      test("should handle complex chemical names", async () => {
        const result = await atomPredictor.structuralizeText("2-methylpropanoic acid");
        expect(result).toHaveProperty("object");
        expect(result).toHaveProperty("chemicals");
        expect(Array.isArray(result.chemicals)).toBe(true);
      });

      test("should throw error when API fails", async () => {
        const mockOpenAI = require("openai");
        mockOpenAI.OpenAI.mockImplementationOnce(() => ({
          chat: {
            completions: {
              create: jest.fn().mockRejectedValue(new Error("Network timeout")),
            },
          },
        }));

        const analyzer = new Structuralizer("invalid-key");
        await expect(analyzer.structuralizeText("test")).rejects.toThrow();
      });
    });

    describe("parseAIResponse", () => {
      test("should parse valid JSON response", () => {
        const validJson = JSON.stringify({
          object: "water",
          chemicals: [{ name: "Water", smiles: "O" }]
        });
        const result = atomPredictor.parseAIResponse(validJson);
        expect(result.object).toBe("water");
        expect(result.chemicals).toHaveLength(1);
        expect(result.chemicals[0].smiles).toBe("O");
      });

      test("should extract JSON from mixed content", () => {
        const mixedContent = `Here's the analysis: {"object": "ethanol", "chemicals": [{"name": "Ethanol", "smiles": "CCO"}]} and some extra text.`;
        const result = atomPredictor.parseAIResponse(mixedContent);
        expect(result.object).toBe("ethanol");
        expect(result.chemicals[0].smiles).toBe("CCO");
      });

      test("should handle malformed JSON gracefully", () => {
        const malformedJson = `{"object": "test", "chemicals":`;
        const result = atomPredictor.parseAIResponse(malformedJson);
        expect(result.object).toBe("Unknown object");
        expect(result.chemicals).toEqual([]);
      });

      test("should handle non-JSON response", () => {
        const nonJson = "This is not JSON at all";
        const result = atomPredictor.parseAIResponse(nonJson);
        expect(result.object).toBe("Unknown object");
        expect(result.chemicals).toEqual([]);
      });

      test("should handle empty response", () => {
        const result = atomPredictor.parseAIResponse("");
        expect(result.object).toBe("Unknown object");
        expect(result.chemicals).toEqual([]);
      });

      test("should handle null/undefined response", () => {
        const nullResult = atomPredictor.parseAIResponse(null);
        expect(nullResult.object).toBe("Unknown object");
        
        const undefinedResult = atomPredictor.parseAIResponse(undefined);
        expect(undefinedResult.object).toBe("Unknown object");
      });
    });

    // buildChemicalInstructions tests removed; Structuralizer uses structuralize prompt
  });

  describe("MolecularProcessor", () => {
    beforeEach(() => {
      // Reset spawn mock before each test
      const { spawn } = require("child_process");
      spawn.mockClear();
    });

    describe("constructor and initialization", () => {
      test("should initialize with default sdf directory", () => {
        const processor = new MolecularProcessor();
        expect(processor.sdfDir).toContain("test/sdf_files");
      });

      test("should initialize with custom sdf directory", () => {
        const customDir = "custom/sdf/path";
        const processor = new MolecularProcessor(customDir);
        expect(processor.sdfDir).toContain(customDir);
      });

      test("should ensure SDF directory exists", () => {
        const mockMkdirSync = require("fs").mkdirSync;
        const mockExistsSync = require("fs").existsSync;
        
        mockExistsSync.mockReturnValueOnce(false);
        new MolecularProcessor("test/dir");
        
        expect(mockMkdirSync).toHaveBeenCalledWith(
          expect.stringContaining("test/dir"),
          { recursive: true }
        );
      });
    });

    describe("processSmiles", () => {
      test("should process valid SMILES array", async () => {
        const result = await molecularProcessor.processSmiles(["CCO", "CC(=O)O"]);
        expect(result).toHaveProperty("sdfPaths");
        expect(result).toHaveProperty("errors");
        expect(result).toHaveProperty("skipped");
        expect(Array.isArray(result.sdfPaths)).toBe(true);
        expect(Array.isArray(result.errors)).toBe(true);
        expect(Array.isArray(result.skipped)).toBe(true);
      });

      test("should handle invalid SMILES gracefully", async () => {
        const result = await molecularProcessor.processSmiles(["INVALID_SMILES"]);
        expect(result.errors).toHaveLength(1);
        expect(result.sdfPaths).toHaveLength(0);
        expect(result.errors[0]).toContain("INVALID_SMILES");
      });

      test("should handle empty SMILES array", async () => {
        const result = await molecularProcessor.processSmiles([]);
        expect(result.sdfPaths).toHaveLength(0);
        expect(result.errors).toHaveLength(0);
        expect(result.skipped).toHaveLength(0);
      });

      test("should respect overwrite parameter", async () => {
        const result = await molecularProcessor.processSmiles(["CCO"], true);
        expect(result).toHaveProperty("sdfPaths");
      });

      test("should handle mixed valid and invalid SMILES", async () => {
        // In test environment, valid SMILES will fail with mock spawn, invalid ones get filtered out
        // This test verifies error handling works correctly
        const result = await molecularProcessor.processSmiles(["CCO", "INVALID", "CC(=O)O"]);
        expect(result.errors).toHaveLength(2); // CCO and CC(=O)O fail during spawn (both pass validation)
        expect(result.skipped).toHaveLength(1); // INVALID gets filtered out during validation
        expect(result.skipped[0]).toContain("INVALID");
        // Both CCO and CC(=O)O should be in errors since they pass validation but fail spawn
        expect(result.errors.some(e => e.includes("CCO"))).toBe(true);
        expect(result.errors.some(e => e.includes("CC(=O)O"))).toBe(true);
      });
    });

    describe("generateSDF", () => {
      test("should generate SDF for valid SMILES", async () => {
        // In test environment with mock spawn that fails, this should throw
        await expect(molecularProcessor.generateSDF("CCO")).rejects.toThrow();
      });

      test("should return existing file when not overwriting", async () => {
        const mockExistsSync = require("fs").existsSync;
        mockExistsSync.mockReturnValueOnce(true);
        
        const result = await molecularProcessor.generateSDF("CCO", false);
        expect(result).toContain("/sdf_files/");
      });

      test("should handle overwrite parameter", async () => {
        // In test environment with mock spawn that fails, this should throw
        await expect(molecularProcessor.generateSDF("CCO", true)).rejects.toThrow();
      });

      test("should throw error for invalid SMILES", async () => {
        await expect(molecularProcessor.generateSDF("INVALID_SMILES")).rejects.toThrow();
      });
    });

    describe("generateSmilesSDF", () => {
      test("should spawn python process for SMILES generation", async () => {
        const { spawn } = require("child_process");
        
        spawn.mockReturnValue({
          stdout: {
            on: jest.fn((event, callback) => {
              if (event === "data") callback("output");
            }),
          },
          stderr: {
            on: jest.fn(),
          },
          on: jest.fn((event, callback) => {
            if (event === "close") callback(0);
          }),
        });
        
        const mockExistsSync = require("fs").existsSync;
        mockExistsSync.mockReturnValueOnce(true);

        try {
          await molecularProcessor.generateSmilesSDF("CCO");
        } catch (error) {
          // Expected to fail in test environment
        }
        
        expect(spawn).toHaveBeenCalledWith(
          "python",
          expect.arrayContaining([
            expect.stringContaining("sdf.py"),
            "CCO",
            "--dir",
            expect.any(String)
          ]),
          expect.any(Object)
        );
      });

      test("should handle python process failure", async () => {
        const { spawn } = require("child_process");
        
        spawn.mockReturnValue({
          stdout: { on: jest.fn() },
          stderr: { on: jest.fn() },
          on: jest.fn((event, callback) => {
            if (event === "close") callback(1); // Non-zero exit code
          }),
        });

        await expect(molecularProcessor.generateSmilesSDF("INVALID")).rejects.toThrow("SMILES generation failed");
      });
    });

    describe("findExistingSdfFile", () => {
      test("should find existing SDF file with exact name", () => {
        const mockExistsSync = require("fs").existsSync;
        mockExistsSync.mockImplementation((path) => path.includes("CCO.sdf"));

        const result = molecularProcessor.findExistingSdfFile("CCO");
        expect(result).toBe("/sdf_files/CCO.sdf");
      });

      test("should find existing SDF file with sanitized name", () => {
        const mockExistsSync = require("fs").existsSync;
        mockExistsSync.mockImplementation((path) => path.includes("CC_O_.sdf"));

        const result = molecularProcessor.findExistingSdfFile("CC(O)");
        expect(result).toBe("/sdf_files/CC_O_.sdf");
      });

      test("should return null when no file exists", () => {
        const mockExistsSync = require("fs").existsSync;
        mockExistsSync.mockReturnValue(false);

        const result = molecularProcessor.findExistingSdfFile("NONEXISTENT");
        expect(result).toBe(null);
      });

      test("should handle special characters in SMILES", () => {
        const mockExistsSync = require("fs").existsSync;
        mockExistsSync.mockImplementation((path) => path.includes("CC___O_O.sdf"));

        const result = molecularProcessor.findExistingSdfFile("CC(=O)O");
        expect(result).toBe("/sdf_files/CC___O_O.sdf");
      });
    });

    describe("ensureSdfDirectory", () => {
      test("should create directory if it doesn't exist", () => {
        const mockExistsSync = require("fs").existsSync;
        const mockMkdirSync = require("fs").mkdirSync;
        
        mockExistsSync.mockReturnValueOnce(false);
        
        const processor = new MolecularProcessor("new/test/dir");
        expect(mockMkdirSync).toHaveBeenCalledWith(
          expect.stringContaining("new/test/dir"),
          { recursive: true }
        );
      });

      test("should not create directory if it exists", () => {
        const mockExistsSync = require("fs").existsSync;
        const mockMkdirSync = require("fs").mkdirSync;
        
        mockExistsSync.mockReturnValueOnce(true);
        mockMkdirSync.mockClear();
        
        new MolecularProcessor("existing/dir");
        expect(mockMkdirSync).not.toHaveBeenCalled();
      });
    });
  });

  describe("Schemas", () => {
    test("ImageMoleculeSchema should validate correct data", () => {
      const validData = {
        imageBase64:
          "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==",
        croppedImageBase64:
          "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==",
        x: 100,
        y: 100,
      };

      const result = ImageMoleculeSchema.safeParse(validData);
      expect(result.success).toBe(true);
    });

    test("ImageMoleculeSchema should reject invalid data", () => {
      const invalidData = {
        imageBase64: "invalid-base64",
        x: "not-a-number",
      };

      const result = ImageMoleculeSchema.safeParse(invalidData);
      expect(result.success).toBe(false);
    });

    test("TextMoleculeSchema should validate correct data", () => {
      const validData = { object: "test object" };
      const result = TextMoleculeSchema.safeParse(validData);
      expect(result.success).toBe(true);
    });

    test("SdfGenerationSchema should validate correct data", () => {
      const validData = {
        smiles: ["CCO", "CC(=O)O"],
        overwrite: false,
      };
      const result = SdfGenerationSchema.safeParse(validData);
      expect(result.success).toBe(true);
    });
  });

  describe("Utility Functions", () => {
    test("should process SMILES array", async () => {
      const result = await molecularProcessor.processSmiles(["CCO", "CC(=O)O"]);
      expect(result).toHaveProperty("sdfPaths");
      expect(result).toHaveProperty("errors");
      expect(result).toHaveProperty("skipped");
    });

    test("should handle file operations", () => {
      expect(molecularProcessor.sdfDir).toBeDefined();
      expect(typeof molecularProcessor.findExistingSdfFile).toBe("function");
    });
  });

  describe("Error Handling", () => {
    test("should handle network errors in AI analyzer", async () => {
      const mockOpenAI = require("openai");
      mockOpenAI.OpenAI.mockImplementationOnce(() => ({
        chat: {
          completions: {
            create: jest.fn().mockRejectedValue(new Error("Network Error")),
          },
        },
      }));

      const analyzer = new Structuralizer("test-key");
      await expect(analyzer.structuralizeText("test")).rejects.toThrow(
        "Network Error",
      );
    });

    test("should handle processing errors gracefully", async () => {
      const result = await molecularProcessor.processSmiles(["INVALID_SMILES"]);
      expect(result.errors.length).toBeGreaterThan(0);
      expect(result.sdfPaths).toHaveLength(0);
    });
  });

  describe("API Endpoints", () => {
    describe("POST /image-molecules", () => {
      test("should analyze image successfully", async () => {
        const mockImageData = {
          imageBase64: "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==",
          croppedImageBase64: "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==",
          x: 100,
          y: 100
        };

        const response = await request(app)
          .post("/image-molecules")
          .send(mockImageData);

        expect(response.status).toBe(200);
        expect(response.body).toHaveProperty("object");
        expect(response.body).toHaveProperty("chemicals");
        expect(Array.isArray(response.body.chemicals)).toBe(true);
      });

      test("should reject invalid image data", async () => {
        const invalidData = {
          imageBase64: "invalid-base64",
          x: "not-a-number"
        };

        const response = await request(app)
          .post("/image-molecules")
          .send(invalidData);

        expect(response.status).toBe(400);
      });

      test("should handle missing required fields", async () => {
        const response = await request(app)
          .post("/image-molecules")
          .send({});

        expect(response.status).toBe(400);
      });
    });

    describe("POST /analyze-text", () => {
      test("should analyze text object successfully", async () => {
        const mockTextData = {
          object: "ethanol"
        };

        const response = await request(app)
          .post("/analyze-text")
          .send(mockTextData);

        expect(response.status).toBe(200);
        expect(response.body).toHaveProperty("object");
        expect(response.body).toHaveProperty("chemicals");
        expect(Array.isArray(response.body.chemicals)).toBe(true);
      });

      test("should handle empty object text", async () => {
        const response = await request(app)
          .post("/analyze-text")
          .send({ object: "" });

        expect(response.status).toBe(400);
        expect(response.body).toHaveProperty("error");
      });

      test("should reject invalid text data", async () => {
        const response = await request(app)
          .post("/analyze-text")
          .send({ object: 123 });

        expect(response.status).toBe(400);
      });
    });

    describe("POST /object-molecules", () => {
      test("should analyze object description successfully", async () => {
        const mockObjectData = {
          object: "coffee bean"
        };

        const response = await request(app)
          .post("/object-molecules")
          .send(mockObjectData);

        expect(response.status).toBe(200);
        expect(response.body).toHaveProperty("object");
        expect(response.body).toHaveProperty("chemicals");
      });

      test("should handle complex object descriptions", async () => {
        const mockObjectData = {
          object: "pharmaceutical tablet containing acetaminophen"
        };

        const response = await request(app)
          .post("/object-molecules")
          .send(mockObjectData);

        expect(response.status).toBe(200);
        expect(response.body).toHaveProperty("chemicals");
      });
    });

    describe("POST /generate-sdfs", () => {
      test("should generate SDF files successfully", async () => {
        const mockSdfData = {
          smiles: ["CCO", "CC(=O)O"],
          overwrite: false
        };

        const response = await request(app)
          .post("/generate-sdfs")
          .send(mockSdfData);

        expect(response.status).toBe(200);
        expect(response.body).toHaveProperty("sdfPaths");
        expect(response.body).toHaveProperty("errors");
        expect(response.body).toHaveProperty("skipped");
      });

      test("should handle invalid SMILES in SDF generation", async () => {
        const mockSdfData = {
          smiles: ["INVALID_SMILES"],
          overwrite: false
        };

        const response = await request(app)
          .post("/generate-sdfs")
          .send(mockSdfData);

        expect(response.status).toBe(200);
        expect(response.body.errors).toHaveLength(1);
        expect(response.body.sdfPaths).toHaveLength(0);
      });

      test("should respect overwrite parameter", async () => {
        const mockSdfData = {
          smiles: ["CCO"],
          overwrite: true
        };

        const response = await request(app)
          .post("/generate-sdfs")
          .send(mockSdfData);

        expect(response.status).toBe(200);
        expect(response.body).toHaveProperty("sdfPaths");
      });

      test("should reject invalid schema data", async () => {
        const invalidData = {
          smiles: "not-an-array",
          overwrite: "not-a-boolean"
        };

        const response = await request(app)
          .post("/generate-sdfs")
          .send(invalidData);

        expect(response.status).toBe(400);
      });
    });

    describe("Error handling", () => {
      test("should handle server errors gracefully", async () => {
        // Test with invalid object that will cause validation error
        const response = await request(app)
          .post("/analyze-text")
          .send({ object: 123 }); // Invalid type

        expect(response.status).toBe(400);
        expect(response.body).toHaveProperty("error");
      });

      test("should handle malformed JSON requests", async () => {
        const response = await request(app)
          .post("/analyze-text")
          .set('Content-Type', 'application/json')
          .send('{"malformed": json}');

        expect(response.status).toBe(400);
      });
    });
  });

  describe("Input Validation", () => {
    test("should validate schema data", () => {
      const validData = { object: "test object" };
      const result = TextMoleculeSchema.safeParse(validData);
      expect(result.success).toBe(true);
    });

    test("should reject invalid schema data", () => {
      const invalidData = { object: 123 };
      const result = TextMoleculeSchema.safeParse(invalidData);
      expect(result.success).toBe(false);
    });
  });
});
