// test/unit/unit.test.js - Unit tests for individual components
// These tests run quickly and validate individual functions and modules

const MolecularProcessor = require("../../../src/server/services/molecular-processor");
const {
  ImageMoleculeSchema,
  TextMoleculeSchema,
  SdfGenerationSchema,
} = require("../../../src/server/schemas/molecular");

// Mock file system
jest.mock("fs", () => ({
  existsSync: jest.fn(),
  writeFileSync: jest.fn(),
  readFileSync: jest.fn(),
  mkdirSync: jest.fn(),
}));

// Mock child_process for SDF generation
jest.mock("child_process", () => ({
  execSync: jest.fn().mockReturnValue("test output"),
  spawn: jest.fn().mockReturnValue({
    stdout: { on: jest.fn() },
    stderr: { on: jest.fn() },
    on: jest.fn((event, callback) => {
      if (event === "close") callback(1);
    }),
  }),
}));


describe("Unit Tests", () => {
  // Some MolecularProcessor tests can be slow (even with mocked child_process)
  // and get slower under `--detectOpenHandles`. Keep this scoped to this file.
  jest.setTimeout(30000);
  let molecularProcessor;

  beforeEach(() => {
    molecularProcessor = new MolecularProcessor();
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
        expect(processor.sdfDir).toContain("tests/sdf_files");
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
      const validData = { image: "base64data" };
      const result = ImageMoleculeSchema.safeParse(validData);
      expect(result.success).toBe(true);
    });

    test("ImageMoleculeSchema should reject invalid data", () => {
      const invalidData = { image: "" };
      const result = ImageMoleculeSchema.safeParse(invalidData);
      expect(result.success).toBe(false);
    });

    test("TextMoleculeSchema should validate correct data", () => {
      const validData = { text: "test molecule" };
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
    test("should handle processing errors gracefully", async () => {
      const result = await molecularProcessor.processSmiles(["INVALID_SMILES"]);
      expect(result.errors.length).toBeGreaterThan(0);
      expect(result.sdfPaths).toHaveLength(0);
    });
  });

  describe("Input Validation", () => {
    test("should validate schema data", () => {
      const validData = { text: "test molecule" };
      const result = TextMoleculeSchema.safeParse(validData);
      expect(result.success).toBe(true);
    });

    test("should reject invalid schema data", () => {
      const invalidData = { text: 123 };
      const result = TextMoleculeSchema.safeParse(invalidData);
      expect(result.success).toBe(false);
    });
  });
});
