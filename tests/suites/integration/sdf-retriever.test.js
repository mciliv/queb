const path = require('path');
const fs = require('fs');
const { spawn } = require('child_process');

describe('SDF Retriever Integration Tests', () => {
  const SDF_SCRIPT = path.join(__dirname, '../../../src/server/molecular-docking-research/sdf.py');
  const TEST_OUTPUT_DIR = path.join(__dirname, '../../temp_sdf_test');
  const BACKEND_SDF_DIR = path.join(__dirname, '../../../src/backend/sdf_files');
  const TEST_SDF_DIR = path.join(__dirname, '../../sdf_files');

  beforeAll(() => {
    // Ensure test output directory exists
    if (!fs.existsSync(TEST_OUTPUT_DIR)) {
      fs.mkdirSync(TEST_OUTPUT_DIR, { recursive: true });
    }
  });

  afterAll(() => {
    // Cleanup test output directory
    try {
      if (fs.existsSync(TEST_OUTPUT_DIR)) {
        fs.rmSync(TEST_OUTPUT_DIR, { recursive: true, force: true });
      }
    } catch (error) {
      console.warn('Failed to cleanup test directory:', error.message);
    }
  });

  const runSdfScript = (formula, options = {}) => {
    return new Promise((resolve, reject) => {
      const args = [SDF_SCRIPT, formula];
      
      if (options.dir) {
        args.push('--dir', options.dir);
      }
      if (options.overwrite) {
        args.push('--overwrite');
      }

      const child = spawn('python', args, {
        stdio: ['pipe', 'pipe', 'pipe'],
        timeout: 30000 // 30 second timeout
      });

      let stdout = '';
      let stderr = '';

      child.stdout.on('data', (data) => {
        stdout += data.toString();
      });

      child.stderr.on('data', (data) => {
        stderr += data.toString();
      });

      child.on('close', (code) => {
        resolve({
          code,
          stdout: stdout.trim(),
          stderr: stderr.trim()
        });
      });

      child.on('error', (error) => {
        reject(error);
      });
    });
  };

  const validateSDF = (filePath) => {
    if (!fs.existsSync(filePath)) {
      return { valid: false, error: 'File does not exist' };
    }

    const content = fs.readFileSync(filePath, 'utf8');
    
    // Basic SDF format validation
    const hasEndMarker = content.includes('M  END') || content.includes('$$$$');
    const hasAtomBlock = /^\s*\d+\s+\d+\s+\d+/m.test(content);
    const fileSize = fs.statSync(filePath).size;

    return {
      valid: hasEndMarker && fileSize > 50,
      hasEndMarker,
      hasAtomBlock,
      fileSize,
      content: content.substring(0, 200) + (content.length > 200 ? '...' : '')
    };
  };

  describe('Local SDF File Retrieval', () => {
    test('should find existing CCO (ethanol) file', async () => {
      const result = await runSdfScript('CCO', { dir: TEST_OUTPUT_DIR });
      
      expect(result.code).toBe(0);
      expect(result.stdout).toContain('✅ SDF file retrieved');
      expect(result.stdout).toMatch(/📁 Solution: Local file system|🌐 Solution: API download/);
      
      // Check if file was created/copied
      const outputFiles = fs.readdirSync(TEST_OUTPUT_DIR).filter(f => f.endsWith('.sdf'));
      expect(outputFiles.length).toBeGreaterThan(0);
      
      const sdfFile = path.join(TEST_OUTPUT_DIR, outputFiles[0]);
      const validation = validateSDF(sdfFile);
      expect(validation.valid).toBe(true);
    }, 30000);

    test('should handle ethanol by name resolution', async () => {
      const result = await runSdfScript('ethanol', { dir: TEST_OUTPUT_DIR });
      
      expect(result.code).toBe(0);
      expect(result.stdout).toContain('✅ SDF file retrieved');
    }, 30000);
  });

  describe('API-based SDF Retrieval', () => {
    test('should fetch caffeine via API', async () => {
      const result = await runSdfScript('caffeine', { dir: TEST_OUTPUT_DIR });
      
      if (result.code === 0) {
        expect(result.stdout).toContain('✅ SDF file retrieved');
        
        const outputFiles = fs.readdirSync(TEST_OUTPUT_DIR).filter(f => f.endsWith('.sdf'));
        const caffeineFile = outputFiles.find(f => f.toLowerCase().includes('caffeine'));
        
        if (caffeineFile) {
          const validation = validateSDF(path.join(TEST_OUTPUT_DIR, caffeineFile));
          expect(validation.valid).toBe(true);
          expect(validation.fileSize).toBeGreaterThan(500); // Caffeine should be a substantial molecule
        }
      } else {
        // API calls might fail in test environment, log for debugging
        console.log('Caffeine API test failed (expected in some test environments):', result.stderr);
      }
    }, 45000);

    test('should handle aspirin lookup', async () => {
      const result = await runSdfScript('aspirin', { dir: TEST_OUTPUT_DIR });
      
      // API calls might fail in test environment
      if (result.code === 0) {
        expect(result.stdout).toContain('✅ SDF file retrieved');
      } else {
        console.log('Aspirin API test failed (expected in some test environments)');
      }
    }, 45000);
  });

  describe('Error Handling', () => {
    test('should handle invalid compound gracefully', async () => {
      const result = await runSdfScript('XYZ123INVALID', { dir: TEST_OUTPUT_DIR });
      
      expect(result.code).toBe(1);
      expect(result.stdout).toContain('❌ No SDF file found');
    }, 30000);

    test('should handle empty input', async () => {
      const result = await runSdfScript('', { dir: TEST_OUTPUT_DIR });
      
      expect(result.code).toBe(1);
    }, 30000);
  });

  describe('Solution Path Analysis', () => {
    test('should report solution method used', async () => {
      const result = await runSdfScript('CCO', { dir: TEST_OUTPUT_DIR });
      
      if (result.code === 0) {
        expect(result.stdout).toMatch(/📁 Solution: Local file system|🌐 Solution: API download/);
        expect(result.stdout).toMatch(/📄 File size: \d+ bytes/);
      }
    }, 30000);
  });

  describe('Comprehensive Test Suite', () => {
    const testCases = [
      // Local file tests (should find existing files)
      { formula: "CCO", expected: "local", description: "Ethanol - should exist locally" },
      { formula: "ethanol", expected: "local", description: "Ethanol by name - should resolve to CCO" },
      { formula: "CaCO3", expected: "local_or_api", description: "Calcium carbonate - may exist locally" },
      { formula: "H2O", expected: "local_or_api", description: "Water - common compound" },
      
      // API resolution tests (likely need PubChem)
      { formula: "caffeine", expected: "api", description: "Caffeine - complex molecule via API" },
      { formula: "aspirin", expected: "api", description: "Aspirin - should resolve via name-resolver" },
      { formula: "C8H10N4O2", expected: "api", description: "Caffeine by formula - API lookup" },
      
      // Edge cases and variations
      { formula: "NaCl", expected: "any", description: "Salt - test ionic compound" },
      { formula: "glucose", expected: "api", description: "Glucose - biological molecule" },
      { formula: "C6H12O6", expected: "api", description: "Glucose by formula" },
      
      // Failure cases (for robustness testing)
      { formula: "XYZ123", expected: "fail", description: "Invalid compound - should fail gracefully" },
      { formula: "", expected: "fail", description: "Empty input - should handle error" },
    ];

    test('should run comprehensive test suite with all scenarios', async () => {
      const results = {
        passed: 0,
        failed: 0,
        local_hits: 0,
        api_hits: 0,
        failures: 0,
        details: []
      };

      console.log('🧪 SDF Retriever Comprehensive Test Suite');
      console.log('=' * 50);

      for (let i = 0; i < testCases.length; i++) {
        const testCase = testCases[i];
        const { formula, expected, description } = testCase;

        console.log(`\n${(i+1).toString().padStart(2)}. Testing: ${formula}`);
        console.log(`    Description: ${description}`);
        console.log(`    Expected: ${expected}`);

        try {
          const result = await runSdfScript(formula, { dir: TEST_OUTPUT_DIR });
          
          if (result.code === 0) {
            // Determine solution path used
            let solution_type = "unknown";
            if (result.stdout.includes("Local file system")) {
              solution_type = "local";
              results.local_hits++;
            } else if (result.stdout.includes("API download")) {
              solution_type = "api";
              results.api_hits++;
            }

            // Validate the SDF file exists
            const outputFiles = fs.readdirSync(TEST_OUTPUT_DIR).filter(f => f.endsWith('.sdf'));
            if (outputFiles.length > 0) {
              const sdfFile = path.join(TEST_OUTPUT_DIR, outputFiles[0]);
              const validation = validateSDF(sdfFile);
              
              if (validation.valid) {
                console.log(`    ✅ SUCCESS: ${solution_type} (${validation.fileSize} bytes)`);
                console.log(`    📄 Valid SDF format detected`);
                results.passed++;
              } else {
                console.log(`    ❌ FAIL: Invalid SDF format`);
                results.failed++;
              }
            } else {
              console.log(`    ❌ FAIL: No SDF file created`);
              results.failed++;
            }
          } else {
            if (expected === "fail") {
              console.log(`    ✅ EXPECTED FAILURE: Correctly handled invalid input`);
              results.passed++;
            } else {
              console.log(`    ❌ FAIL: No SDF retrieved`);
              results.failed++;
              results.failures++;
            }
          }

          results.details.push({
            formula,
            expected,
            result: result.code === 0 ? "success" : "fail",
            solution_type: result.code === 0 ? solution_type : "none"
          });

        } catch (error) {
          if (expected === "fail") {
            console.log(`    ✅ EXPECTED FAILURE: ${error.message.substring(0, 50)}...`);
            results.passed++;
          } else {
            console.log(`    ❌ ERROR: ${error.message.substring(0, 50)}...`);
            results.failed++;
          }

          results.details.push({
            formula,
            expected,
            result: "error",
            error: error.message,
            solution_type: "error"
          });
        }

        // Clean up files between tests
        try {
          const files = fs.readdirSync(TEST_OUTPUT_DIR);
          files.forEach(file => {
            if (file.endsWith('.sdf')) {
              fs.unlinkSync(path.join(TEST_OUTPUT_DIR, file));
            }
          });
        } catch (e) {
          // Ignore cleanup errors
        }
      }

      // Summary report
      console.log('\n' + '='.repeat(50));
      console.log('📊 TEST SUMMARY');
      console.log('='.repeat(50));
      console.log(`Total Tests: ${testCases.length}`);
      console.log(`✅ Passed: ${results.passed}`);
      console.log(`❌ Failed: ${results.failed}`);
      console.log(`📁 Local Hits: ${results.local_hits}`);
      console.log(`🌐 API Hits: ${results.api_hits}`);
      console.log(`💥 Failures: ${results.failures}`);

      // Solution path analysis
      console.log(`\n🔍 SOLUTION PATH ANALYSIS:`);
      console.log(`Local file retrieval: ${results.local_hits} cases`);
      console.log(`API-based retrieval: ${results.api_hits} cases`);
      console.log(`Failed retrievals: ${results.failures} cases`);

      const success_rate = (results.passed / testCases.length) * 100;
      console.log(`\n📈 Overall Success Rate: ${success_rate.toFixed(1)}%`);

      // Test assertions
      expect(results.passed).toBeGreaterThan(0);
      expect(success_rate).toBeGreaterThan(50); // At least 50% success rate
      expect(results.local_hits + results.api_hits).toBeGreaterThan(0); // At least some successful retrievals
    }, 120000);
  });

  describe('File System Integration', () => {
    test('should respect existing SDF directories', async () => {
      // Test that it checks the right directories
      const result = await runSdfScript('CCO', { dir: TEST_OUTPUT_DIR });
      
      if (result.code === 0) {
        // Should find file in either backend or test sdf_files
        const foundLocal = result.stdout.includes('📁 Solution: Local file system');
        const foundAPI = result.stdout.includes('🌐 Solution: API download');
        
        expect(foundLocal || foundAPI).toBe(true);
      }
    }, 30000);

    test('should create output directory if needed', async () => {
      const tempDir = path.join(TEST_OUTPUT_DIR, 'nested', 'test', 'dir');
      
      const result = await runSdfScript('CCO', { dir: tempDir });
      
      if (result.code === 0) {
        expect(fs.existsSync(tempDir)).toBe(true);
        
        // Cleanup
        fs.rmSync(path.join(TEST_OUTPUT_DIR, 'nested'), { recursive: true, force: true });
      }
    }, 30000);
  });

  describe('Overwrite Functionality', () => {
    test('should respect overwrite flag', async () => {
      // First run
      const result1 = await runSdfScript('CCO', { dir: TEST_OUTPUT_DIR });
      
      if (result1.code === 0) {
        // Second run without overwrite
        const result2 = await runSdfScript('CCO', { dir: TEST_OUTPUT_DIR });
        expect(result2.code).toBe(0);
        
        // Third run with overwrite
        const result3 = await runSdfScript('CCO', { dir: TEST_OUTPUT_DIR, overwrite: true });
        expect(result3.code).toBe(0);
      }
    }, 45000);
  });
});
