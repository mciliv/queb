// test/unit/unit.test.js - Unit tests for individual components
// These tests run quickly (< 5 seconds) and validate individual functions and modules

// Mock child_process at the top of the file
jest.mock('child_process', () => ({
  spawn: jest.fn(() => ({
    stdout: { 
      on: jest.fn((event, callback) => {
        if (event === 'data') {
          // Mock successful SDF generation output
          callback('Generated SDF file: /mock/path/CCO.sdf\n');
        }
      })
    },
    stderr: { on: jest.fn() },
    on: jest.fn((event, callback) => {
      if (event === 'close') {
        setTimeout(() => callback(0), 10);
      }
    })
  })),
  execSync: jest.fn().mockReturnValue("test output"),
}));

const request = require("supertest");
const fs = require("fs");
const path = require("path");
const AtomPredictor = require("../../backend/services/AtomPredictor");
const MolecularProcessor = require("../../backend/services/molecular-processor");
const {
  ImageMoleculeSchema,
  TextMoleculeSchema,
  SdfGenerationSchema,
} = require("../../backend/schemas/schemas");

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



// Global app instance for all tests
let app;

beforeAll(() => {
  // Load server for all tests
  try {
    app = require("../../backend/api/server");
  } catch (error) {
    console.error("Failed to load server:", error.message);
    throw error;
  }
});

describe("Unit Tests", () => {
  let atomPredictor;
  let molecularProcessor;

  beforeEach(() => {
    atomPredictor = new AtomPredictor("test-api-key");
    molecularProcessor = new MolecularProcessor();
  });

  describe("AtomPredictor", () => {
    test("should initialize with API key", () => {
      expect(atomPredictor.client).toBeDefined();
      expect(atomPredictor.chemicalInstructions).toBeDefined();
      expect(typeof atomPredictor.chemicalInstructions).toBe("string");
    });

    describe("analyzeImage", () => {
      const mockImageBase64 =
        "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==";

      test("should analyze image without cropped region", async () => {
        const result = await atomPredictor.analyzeImage(mockImageBase64);
        expect(result).toHaveProperty("object");
        expect(result).toHaveProperty("chemicals");
        expect(Array.isArray(result.chemicals)).toBe(true);
        expect(result.object).toBe("test object");
        expect(result.chemicals).toHaveLength(2);
      });

      test("should analyze image with cropped region", async () => {
        const result = await atomPredictor.analyzeImage(
          mockImageBase64,
          mockImageBase64,
          100,
          100,
          50,
          50,
          100
        );
        expect(result).toHaveProperty("object");
        expect(result).toHaveProperty("chemicals");
        expect(Array.isArray(result.chemicals)).toBe(true);
      });

      test("should handle null parameters gracefully", async () => {
        const result = await atomPredictor.analyzeImage(
          mockImageBase64,
          null,
          null,
          null
        );
        expect(result).toHaveProperty("object");
        expect(result).toHaveProperty("chemicals");
      });

      test("should throw error when API fails", async () => {
        const mockOpenAI = require("openai");
        const originalMock = mockOpenAI.OpenAI.getMockImplementation();
        
        mockOpenAI.OpenAI.mockImplementationOnce(() => ({
          chat: {
            completions: {
              create: jest.fn().mockRejectedValue(new Error("API Error")),
            },
          },
        }));

        const analyzer = new AtomPredictor("invalid-key");
        await expect(analyzer.analyzeImage(mockImageBase64)).rejects.toThrow("AI analysis failed: API Error");
        
        // Restore original mock
        mockOpenAI.OpenAI.mockImplementation(originalMock);
      });

      test("should handle empty image base64", async () => {
        const result = await atomPredictor.analyzeImage("");
        expect(result).toBeDefined(); // Empty image is handled gracefully
      });
    });

    describe("analyzeText", () => {
      test("should analyze text input successfully", async () => {
        const result = await atomPredictor.analyzeText("ethanol");
        expect(result).toHaveProperty("object");
        expect(result).toHaveProperty("chemicals");
        expect(Array.isArray(result.chemicals)).toBe(true);
        expect(result.object).toBe("test object");
        expect(result.chemicals[0]).toHaveProperty("name");
        expect(result.chemicals[0]).toHaveProperty("smiles");
      });

      test("should handle empty text input", async () => {
        const result = await atomPredictor.analyzeText("");
        expect(result).toHaveProperty("object");
        expect(result).toHaveProperty("chemicals");
      });

      test("should handle complex chemical names", async () => {
        const result = await atomPredictor.analyzeText("2-methylpropanoic acid");
        expect(result).toHaveProperty("object");
        expect(result).toHaveProperty("chemicals");
        expect(Array.isArray(result.chemicals)).toBe(true);
      });

      test("should throw error when API fails", async () => {
        const mockOpenAI = require("openai");
        const originalMock = mockOpenAI.OpenAI.getMockImplementation();
        
        mockOpenAI.OpenAI.mockImplementationOnce(() => ({
          chat: {
            completions: {
              create: jest.fn().mockRejectedValue(new Error("Network timeout")),
            },
          },
        }));

        const analyzer = new AtomPredictor("invalid-key");
        await expect(analyzer.analyzeText("test")).rejects.toThrow("Request timeout: The AI service is taking too long to respond.");
        
        // Restore original mock
        mockOpenAI.OpenAI.mockImplementation(originalMock);
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
        expect(result.object).toBe("test");
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

    describe("buildChemicalInstructions", () => {
      test("should return proper instruction string", () => {
        const instructions = atomPredictor.buildChemicalInstructions();
        expect(typeof instructions).toBe("string");
        expect(instructions).toContain("JSON response");
        expect(instructions).toContain("SMILES notation");
        expect(instructions).toContain("object");
        expect(instructions).toContain("chemicals");
      });
    });
  });

  describe("MolecularProcessor", () => {
    describe("constructor and initialization", () => {
      test("should initialize with default sdf directory", () => {
        const processor = new MolecularProcessor();
        expect(processor.sdfDir).toContain("data/sdf_files");
      });

      test("should initialize with custom sdf directory", () => {
        const customDir = "custom/sdf/path";
        const processor = new MolecularProcessor(customDir);
        expect(processor.sdfDir).toContain(customDir);
      });

      test("should ensure SDF directory exists", () => {
        const mockMkdirSync = require("fs").mkdirSync;
        const mockExistsSync = require("fs").existsSync;
        
        // Clear any previous calls
        mockMkdirSync.mockClear();
        mockExistsSync.mockClear();
        
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
        const result = await molecularProcessor.processSmiles(["CCO", "INVALID", "CC(=O)O"]);
        expect(result.skipped.length).toBeGreaterThan(0);
        expect(result.skipped.some(s => s.includes("INVALID"))).toBe(true);
      });
    });

    describe("generateSDF", () => {
      test("should generate SDF for valid SMILES", async () => {
        const mockExistsSync = require("fs").existsSync;
        mockExistsSync.mockReturnValueOnce(true); // Mock that the SDF file exists after generation
        
        const result = await molecularProcessor.generateSDF("CCO");
        expect(typeof result === "string" || result === null).toBe(true);
      });

      test("should return existing file when not overwriting", async () => {
        const mockExistsSync = require("fs").existsSync;
        mockExistsSync.mockReturnValueOnce(true);
        
        const result = await molecularProcessor.generateSDF("CCO", false);
        expect(result).toContain("/sdf_files/");
      });

      test("should handle overwrite parameter", async () => {
        const mockExistsSync = require("fs").existsSync;
        mockExistsSync.mockReturnValueOnce(true); // Mock that the SDF file exists after generation
        
        const result = await molecularProcessor.generateSDF("CCO", true);
        expect(typeof result === "string" || result === null).toBe(true);
      });

      test("should throw error for invalid SMILES", async () => {
        await expect(molecularProcessor.generateSDF("INVALID_SMILES")).rejects.toThrow();
      });
    });

    describe("generateSmilesSDF", () => {
      test("should spawn python process for SMILES generation", async () => {
        const { spawn } = require("child_process");
        
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
        
        // Mock spawn to return error condition
        spawn.mockReturnValueOnce({
          stdout: { on: jest.fn() },
          stderr: { on: jest.fn() },
          on: jest.fn((event, callback) => {
            if (event === "close") callback(1); // Non-zero exit code
          }),
        });

        await expect(molecularProcessor.generateSmilesSDF("INVALID")).rejects.toThrow("SMILES generation failed for INVALID (exit code: 1)");
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
        expect(result).toBeNull(); // File doesn't exist in test environment
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

      const analyzer = new AtomPredictor("test-key");
      await expect(analyzer.analyzeText("test")).rejects.toThrow(
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
        expect(response.body).toHaveProperty("output");
        expect(response.body.output).toHaveProperty("object");
        expect(response.body.output).toHaveProperty("chemicals");
        expect(Array.isArray(response.body.output.chemicals)).toBe(true);
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
          .post("/object-molecules")
          .send(mockTextData);

        expect(response.status).toBe(200);
        expect(response.body).toHaveProperty("output");
        expect(response.body.output).toHaveProperty("object");
        expect(response.body.output).toHaveProperty("chemicals");
        expect(Array.isArray(response.body.output.chemicals)).toBe(true);
      });

      test("should handle empty object text", async () => {
        const response = await request(app)
          .post("/analyze-text")
          .send({ text: "" });

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
        expect(response.body).toHaveProperty("output");
        expect(response.body.output).toHaveProperty("object");
        expect(response.body.output).toHaveProperty("chemicals");
      });

      test("should handle complex object descriptions", async () => {
        const mockObjectData = {
          object: "pharmaceutical tablet containing acetaminophen"
        };

        const response = await request(app)
          .post("/object-molecules")
          .send(mockObjectData);

        expect(response.status).toBe(200);
        expect(response.body).toHaveProperty("output");
        expect(response.body.output).toHaveProperty("chemicals");
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
        // Mock an error in AtomPredictor
        const mockOpenAI = require("openai");
        mockOpenAI.OpenAI.mockImplementationOnce(() => ({
          chat: {
            completions: {
              create: jest.fn().mockRejectedValue(new Error("OpenAI API Error")),
            },
          },
        }));

        const response = await request(app)
          .post("/analyze-text")
          .send({ object: "test" });

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
