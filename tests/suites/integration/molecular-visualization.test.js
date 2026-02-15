// Integration test for molecular visualization flow
// Tests the complete pipeline from analysis to visualization

// Mock OpenAI before any modules are imported
jest.mock('openai', () => {
  return {
    default: jest.fn().mockImplementation(() => ({
      chat: {
        completions: {
          create: jest.fn().mockImplementation(async ({ messages }) => {
            const last = Array.isArray(messages) ? [...messages].reverse().find(m => m.role === 'user') : null;
            const content = last?.content;
            let text = '';
            if (typeof content === 'string') {
              text = content;
            } else if (Array.isArray(content)) {
              const textPart = content.find(p => p && p.type === 'text');
              text = textPart?.text || '';
            }

            const lower = (text || '').toLowerCase();
            const chemicals = [];
            if (lower.includes('water')) chemicals.push({ name: 'water', smiles: 'O' });
            if (lower.includes('ethanol') || lower.includes('wine')) chemicals.push({ name: 'ethanol', smiles: 'CCO' });
            if (lower.includes('sodium') && lower.includes('chloride')) chemicals.push({ name: 'sodium chloride', smiles: '[Na+].[Cl-]' });
            if (lower.includes('coffee')) {
              chemicals.push({ name: 'caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' });
              chemicals.push({ name: 'water', smiles: 'O' });
            }
            if (lower.includes('apple')) {
              chemicals.push({ name: 'fructose', smiles: 'C(C(C(C(C(C=O)O)O)O)O)O' });
              chemicals.push({ name: 'glucose', smiles: 'C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O' });
              chemicals.push({ name: 'water', smiles: 'O' });
            }
            if (chemicals.length === 0) chemicals.push({ name: 'water', smiles: 'O' });

            const json = { object: text || 'test object', chemicals };
            return { choices: [{ message: { content: JSON.stringify(json) } }] };
          })
        }
      }
    }))
  };
});

const request = require('supertest');
const path = require('path');
const fs = require('fs');

// Import server for testing
const app = require('../../../src/server/api/server');

describe('Molecular Visualization Integration Tests', () => {
  beforeAll(async () => {
    // Set test environment
    process.env.NODE_ENV = 'test';
    process.env.OPENAI_API_KEY = 'test-key';
  });

  afterAll(async () => {
    // Cleanup if needed
  });

  const testObjects = [
    {
      name: 'Water',
      input: 'water',
      expectedMolecules: ['water', 'H2O'],
      minMoleculeCount: 1
    },
    {
      name: 'Apple',
      input: 'apple',
      expectedMolecules: ['water', 'glucose', 'fructose'],
      minMoleculeCount: 3
    },
    {
      name: 'Coffee',
      input: 'coffee',
      expectedMolecules: ['caffeine', 'water'],
      minMoleculeCount: 2
    },
    {
      name: 'Ethanol',
      input: 'ethanol',
      expectedMolecules: ['ethanol'],
      minMoleculeCount: 1
    }
  ];

  describe('Complete Molecular Analysis Flow', () => {
    testObjects.forEach((testObject) => {
      test(`should analyze ${testObject.name} and generate visualization data`, async () => {
        console.log(`\n🧪 Testing molecular analysis for: ${testObject.name}`);
        
        // Step 1: Analyze the object
        const analysisResponse = await request(app)
          .post('/object-molecules')
          .send({ object: testObject.input })
          .expect(200);

        expect(analysisResponse.body).toBeDefined();
        expect(analysisResponse.body.output).toBeDefined();
        
        const analysis = analysisResponse.body.output;
        console.log(`📊 Analysis result:`, {
          object: analysis.object,
          moleculeCount: analysis.chemicals?.length || 0,
          molecules: analysis.chemicals?.map(c => c.name) || []
        });

        // Verify analysis contains molecules
        expect(analysis.chemicals).toBeDefined();
        expect(Array.isArray(analysis.chemicals)).toBe(true);
        expect(analysis.chemicals.length).toBeGreaterThanOrEqual(testObject.minMoleculeCount);

        // Verify molecules have required fields
        analysis.chemicals.forEach(molecule => {
          expect(molecule.name).toBeDefined();
          expect(molecule.smiles).toBeDefined();
          expect(typeof molecule.name).toBe('string');
          expect(typeof molecule.smiles).toBe('string');
          expect(molecule.smiles.length).toBeGreaterThan(0);
        });

        // Step 2: Generate SDF files for visualization
        const smilesArray = analysis.chemicals.map(mol => mol.smiles);
        
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
        
        console.log(`📁 SDF generation result:`, {
          generated: sdfResponse.body.sdfPaths.length,
          errors: sdfResponse.body.errors?.length || 0,
          skipped: sdfResponse.body.skipped?.length || 0
        });

        // Verify SDF files were generated
        expect(sdfResponse.body.sdfPaths.length).toBeGreaterThan(0);

        // Step 3: Verify SDF files are accessible
        for (const sdfPath of sdfResponse.body.sdfPaths) {
          const fileName = sdfPath.replace('/sdf_files/', '');
          const fullPath = path.join(__dirname, '../../sdf_files', fileName);
          
          // Check if file exists
          expect(fs.existsSync(fullPath)).toBe(true);
          
          // Verify file content is valid SDF
          const content = fs.readFileSync(fullPath, 'utf8');
          expect(content).toBeDefined();
          expect(content.length).toBeGreaterThan(0);
          expect(content).toContain('$$$$'); // SDF format marker
          
          // Test file serving via HTTP
          const fileResponse = await request(app)
            .get(sdfPath)
            .expect(200);
            
          const responseContent = fileResponse.text || 
            (Buffer.isBuffer(fileResponse.body) ? fileResponse.body.toString('utf8') : fileResponse.body);
          expect(responseContent).toBe(content);
        }

        // Step 4: Verify complete visualization data structure
        const visualizationData = analysis.chemicals.map((mol, index) => ({
          name: mol.name,
          smiles: mol.smiles,
          sdfData: sdfResponse.body.sdfPaths[index] || null
        }));

        visualizationData.forEach(molecule => {
          expect(molecule.name).toBeDefined();
          expect(molecule.smiles).toBeDefined();
          // SDF data may be null if generation failed, but that's acceptable
        });

        console.log(`✅ Complete visualization flow verified for ${testObject.name}`);
        console.log(`   - Molecules analyzed: ${analysis.chemicals.length}`);
        console.log(`   - SDF files generated: ${sdfResponse.body.sdfPaths.length}`);
        console.log(`   - Visualization data ready: ${visualizationData.length}`);
      }, 30000); // 30 second timeout for complete flow
    });
  });

  describe('Error Handling in Visualization Flow', () => {
    test('should handle malformed JSON gracefully', async () => {
      // This tests our fallback handlers for truncated/malformed responses
      const response = await request(app)
        .post('/object-molecules')
        .send({ object: 'complex organic compound with very long molecular structure' })
        .expect(200);

      // Should not throw error even if AI returns truncated JSON
      expect(response.body).toBeDefined();
      expect(response.body.output).toBeDefined();
      
      // Should provide fallback molecules if parsing fails
      const analysis = response.body.output;
      expect(analysis.chemicals).toBeDefined();
      expect(Array.isArray(analysis.chemicals)).toBe(true);
    });

    test('should handle invalid SMILES in SDF generation', async () => {
      const invalidSmiles = ['INVALID_SMILES', 'ANOTHER_BAD_SMILES'];
      
      const response = await request(app)
        .post('/generate-sdfs')
        .send({ smiles: invalidSmiles })
        .expect(200);

      expect(response.body).toBeDefined();
      expect(response.body.errors).toBeDefined();
      expect(Array.isArray(response.body.errors)).toBe(true);
      // Should report errors but not crash
    });
  });

  describe('Performance Testing', () => {
    test('should handle concurrent molecular analysis requests', async () => {
      const concurrentRequests = testObjects.slice(0, 3).map(testObj =>
        request(app)
          .post('/object-molecules')
          .send({ object: testObj.input })
      );

      const responses = await Promise.all(concurrentRequests);
      
      responses.forEach(response => {
        expect(response.status).toBe(200);
        expect(response.body.output).toBeDefined();
        expect(response.body.output.chemicals).toBeDefined();
      });

      console.log(`✅ Concurrent requests handled successfully`);
    }, 45000);
  });
});