// test/integration/full-pipeline-visualization.test.js 
// Complete pipeline test: LLM Analysis → SDF Generation → Molecular Visualization
// Uses the same chemical test data as existing LLM endpoint tests

const request = require('supertest');
const fs = require('fs');
const path = require('path');

// Import test data from existing molecular accuracy tests
const REFERENCE_MATERIALS = {
  // Basic chemicals - core accuracy test
  basics: {
    "water": {
      expectedSMILES: ["O"],
      requiredNames: ["water"],
      forbiddenSMILES: ["H2O", "HOH"],
      accuracy_threshold: 95
    },
    "ethanol": {
      expectedSMILES: ["CCO"],
      requiredNames: ["ethanol"],
      forbiddenSMILES: ["C2H6O", "C2H5OH"],
      accuracy_threshold: 95
    },
    "sodium chloride": {
      expectedSMILES: ["[Na+].[Cl-]", "[Na+]", "[Cl-]"],
      requiredNames: ["sodium", "chloride"],
      forbiddenSMILES: ["NaCl"],
      accuracy_threshold: 85
    }
  },

  // Common beverages - realistic composition testing
  beverages: {
    "red wine": {
      expectedSMILES: ["CCO", "O", "OC(C(O)C(O)=O)C(O)=O"],
      requiredNames: ["ethanol", "water", "tartaric"],
      mayContain: ["glucose", "fructose", "malic", "tannin"],
      minComponents: 3,
      accuracy_threshold: 75
    },
    "black coffee": {
      expectedSMILES: ["O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"],
      requiredNames: ["water", "caffeine"],
      mayContain: ["chlorogenic", "acid"],
      minComponents: 2,
      accuracy_threshold: 80
    }
  },

  // Biological materials - challenging but realistic
  biological: {
    "fresh apple": {
      expectedSMILES: ["O", "C(C(C(C(C(C=O)O)O)O)O)O"],
      requiredNames: ["water", "fructose", "glucose"],
      mayContain: ["malic", "cellulose", "pectin"],
      minComponents: 4,
      accuracy_threshold: 65
    }
  }
};

describe('Full Pipeline Molecular Visualization Tests', () => {
  let app;
  let hasApiKey;

  beforeAll(async () => {
    // Import server for testing
    app = require('../../backend/api/server');
    process.env.NODE_ENV = 'test';
    
    // Check if we have OpenAI API key for testing
    hasApiKey = !!process.env.OPENAI_API_KEY;
    
    if (!hasApiKey) {
      console.log("⚠️  Skipping full pipeline tests - no OPENAI_API_KEY set");
      console.log("   Set OPENAI_API_KEY environment variable to run complete pipeline tests");
    }
  });

  afterAll(async () => {
    // Cleanup if needed
  });

  describe('Complete Pipeline Flow: Analysis → SDF → Visualization', () => {
    // Test each category of reference materials
    Object.entries(REFERENCE_MATERIALS).forEach(([category, materials]) => {
      describe(`${category.charAt(0).toUpperCase() + category.slice(1)} Materials`, () => {
        Object.entries(materials).forEach(([materialName, testData]) => {
          test(`should complete full pipeline for ${materialName}`, async () => {
            if (!hasApiKey) {
              console.log(`⚠️  Skipping ${materialName} - no API key`);
              return;
            }
            console.log(`\n🧪 Testing full pipeline for: ${materialName}`);
            console.log(`   Category: ${category}`);
            console.log(`   Expected components: ${testData.minComponents || testData.expectedSMILES.length}`);

            // STEP 1: AI Analysis - Object → Molecules
            console.log('   📊 Step 1: AI Analysis...');
            const analysisResponse = await request(app)
              .post('/object-molecules')
              .send({ object: materialName })
              .expect(200);

            expect(analysisResponse.body).toBeDefined();
            expect(analysisResponse.body.output).toBeDefined();
            
            const analysis = analysisResponse.body.output;
            const molecules = analysis.chemicals || [];

            console.log(`   ✅ Analysis complete: ${molecules.length} molecules found`);
            molecules.forEach((mol, i) => {
              console.log(`      ${i + 1}. ${mol.name} (${mol.smiles?.substring(0, 30)}${mol.smiles?.length > 30 ? '...' : ''})`);
            });

            // Validate molecular analysis quality
            expect(molecules).toBeDefined();
            expect(Array.isArray(molecules)).toBe(true);
            expect(molecules.length).toBeGreaterThan(0);

            // Check minimum components if specified
            if (testData.minComponents) {
              expect(molecules.length).toBeGreaterThanOrEqual(testData.minComponents);
            }

            // Validate molecule structure
            molecules.forEach(molecule => {
              expect(molecule.name).toBeDefined();
              expect(molecule.smiles).toBeDefined();
              expect(typeof molecule.name).toBe('string');
              expect(typeof molecule.smiles).toBe('string');
              expect(molecule.smiles.length).toBeGreaterThan(0);
            });

            // Check for forbidden SMILES patterns
            if (testData.forbiddenSMILES) {
              molecules.forEach(mol => {
                testData.forbiddenSMILES.forEach(forbidden => {
                  expect(mol.smiles).not.toBe(forbidden);
                });
              });
            }

            // STEP 2: SDF Generation - SMILES → 3D Structure Files
            console.log('   📁 Step 2: SDF Generation...');
            const smilesArray = molecules.map(mol => mol.smiles).filter(Boolean);
            
            const sdfResponse = await request(app)
              .post('/generate-sdfs')
              .send({ 
                smiles: smilesArray,
                overwrite: false 
              })
              .expect(200);

            expect(sdfResponse.body).toBeDefined();
            expect(sdfResponse.body.sdfPaths).toBeDefined();
            expect(Array.isArray(sdfResponse.body.sdfPaths)).toBe(true);
            
            console.log(`   ✅ SDF generation: ${sdfResponse.body.sdfPaths.length} files created`);
            if (sdfResponse.body.errors && sdfResponse.body.errors.length > 0) {
              console.log(`   ⚠️  SDF errors: ${sdfResponse.body.errors.length}`);
              sdfResponse.body.errors.slice(0, 3).forEach(error => {
                console.log(`      - ${error.substring(0, 80)}...`);
              });
            }

            // Validate SDF generation
            expect(sdfResponse.body.sdfPaths.length).toBeGreaterThan(0);

            // STEP 3: File Accessibility - Verify files can be served
            console.log('   🔗 Step 3: File Accessibility...');
            let accessibleFiles = 0;
            const fileAccessResults = [];

            for (const sdfPath of sdfResponse.body.sdfPaths.slice(0, 5)) { // Test first 5 files
              try {
                const fileResponse = await request(app)
                  .get(sdfPath)
                  .expect(200);
                
                expect(fileResponse.text).toBeDefined();
                expect(fileResponse.text.length).toBeGreaterThan(0);
                expect(fileResponse.text).toContain('$$$$'); // SDF format marker
                
                accessibleFiles++;
                fileAccessResults.push({ path: sdfPath, accessible: true });
              } catch (error) {
                fileAccessResults.push({ path: sdfPath, accessible: false, error: error.message });
              }
            }

            console.log(`   ✅ File access: ${accessibleFiles}/${Math.min(sdfResponse.body.sdfPaths.length, 5)} files accessible`);
            
            // Require at least 70% of files to be accessible
            const accessibilityRate = accessibleFiles / Math.min(sdfResponse.body.sdfPaths.length, 5);
            expect(accessibilityRate).toBeGreaterThanOrEqual(0.7);

            // STEP 4: Visualization Data Structure
            console.log('   🎯 Step 4: Visualization Data...');
            const visualizationData = molecules.map((mol, index) => ({
              name: mol.name,
              smiles: mol.smiles,
              sdfData: sdfResponse.body.sdfPaths[index] ? `file://${sdfResponse.body.sdfPaths[index]}` : null
            }));

            // Validate complete visualization data
            expect(visualizationData.length).toBe(molecules.length);
            visualizationData.forEach(item => {
              expect(item.name).toBeDefined();
              expect(item.smiles).toBeDefined();
              // sdfData may be null if SDF generation failed, which is acceptable
            });

            console.log(`   ✅ Visualization data ready for ${visualizationData.length} molecules`);

            // STEP 5: Quality Assessment
            console.log('   📈 Step 5: Quality Assessment...');
            
            // Check for required chemical names (if specified)
            if (testData.requiredNames) {
              const foundNames = molecules.map(mol => mol.name.toLowerCase());
              const requiredFound = testData.requiredNames.filter(required => 
                foundNames.some(found => found.includes(required.toLowerCase()))
              );
              
              console.log(`   📋 Required chemicals found: ${requiredFound.length}/${testData.requiredNames.length}`);
              console.log(`      Required: ${testData.requiredNames.join(', ')}`);
              console.log(`      Found: ${requiredFound.join(', ')}`);
              
              // For basic materials, require high accuracy
              if (category === 'basics') {
                expect(requiredFound.length).toBeGreaterThanOrEqual(Math.ceil(testData.requiredNames.length * 0.8));
              }
            }

            // Check SMILES quality
            const validSmilesCount = molecules.filter(mol => {
              return mol.smiles && 
                     mol.smiles.length > 0 && 
                     mol.smiles.length < 300 && // Reasonable length
                     !/^[A-Z][a-z]?\d*([A-Z][a-z]?\d*)*$/.test(mol.smiles) || // Not chemical formula
                     ["O", "C", "N", "S", "P", "H"].includes(mol.smiles); // Except simple elements
            }).length;

            const smilesQuality = (validSmilesCount / molecules.length) * 100;
            console.log(`   🧬 SMILES quality: ${smilesQuality.toFixed(1)}% (${validSmilesCount}/${molecules.length})`);
            
            // Require reasonable SMILES quality
            expect(smilesQuality).toBeGreaterThanOrEqual(70);

            // Final Pipeline Success Report
            console.log(`   🎉 Pipeline completed successfully for ${materialName}:`);
            console.log(`      - Molecules analyzed: ${molecules.length}`);
            console.log(`      - SDF files generated: ${sdfResponse.body.sdfPaths.length}`);
            console.log(`      - Files accessible: ${accessibleFiles}`);
            console.log(`      - Visualization ready: ${visualizationData.length} items`);
            console.log(`      - SMILES quality: ${smilesQuality.toFixed(1)}%`);

          }, 45000); // 45 second timeout for complete pipeline
        });
      });
    });
  });

  describe('Pipeline Error Handling', () => {
    test('should handle invalid object gracefully', async () => {
      if (!hasApiKey) {
        console.log('⚠️  Skipping error handling test - no API key');
        return;
      }
      console.log('\n🧪 Testing error handling with invalid input...');
      
      const response = await request(app)
        .post('/object-molecules')
        .send({ object: 'completely-invalid-nonsense-object-xyz123' })
        .expect(200);

      // Should not crash, should provide fallback data
      expect(response.body.output).toBeDefined();
      expect(response.body.output.chemicals).toBeDefined();
      expect(Array.isArray(response.body.output.chemicals)).toBe(true);
      
      console.log(`   ✅ Error handling: Returned ${response.body.output.chemicals.length} fallback molecules`);
    });

    test('should handle SDF generation failures gracefully', async () => {
      console.log('\n🧪 Testing SDF generation error handling...');
      
      const invalidSmiles = ['INVALID_SMILES_123', 'ANOTHER_BAD_SMILES'];
      
      const response = await request(app)
        .post('/generate-sdfs')
        .send({ smiles: invalidSmiles })
        .expect(200);

      expect(response.body).toBeDefined();
      expect(response.body.errors).toBeDefined();
      expect(Array.isArray(response.body.errors)).toBe(true);
      
      console.log(`   ✅ SDF error handling: ${response.body.errors.length} errors reported gracefully`);
    });
  });

  describe('Pipeline Performance', () => {
    test('should complete pipeline within reasonable time', async () => {
      if (!hasApiKey) {
        console.log('⚠️  Skipping performance test - no API key');
        return;
      }
      console.log('\n🧪 Testing pipeline performance...');
      
      const startTime = Date.now();
      
      // Test with simple object for speed
      const analysisResponse = await request(app)
        .post('/object-molecules')
        .send({ object: 'water' })
        .expect(200);
      
      const molecules = analysisResponse.body.output?.chemicals || [];
      
      if (molecules.length > 0) {
        const smilesArray = molecules.map(mol => mol.smiles).filter(Boolean);
        
        await request(app)
          .post('/generate-sdfs')
          .send({ smiles: smilesArray })
          .expect(200);
      }
      
      const totalTime = Date.now() - startTime;
      console.log(`   ⏱️  Pipeline completed in ${totalTime}ms`);
      
      // Should complete within 30 seconds
      expect(totalTime).toBeLessThan(30000);
    });
  });
});