/**
 * Complex Chemical Analysis Integration Test
 *
 * Tests the full chemical analysis pipeline using complex food matrices as test cases.
 * Parameterized tests validate the same functionality across multiple complex chemical inputs
 * (cacao, coffee, wine, etc.) to ensure robustness with real-world chemical compositions.
 *
 * This test validates:
 * - Full Express app pipeline with mocked AI service
 * - Chemical analysis of complex food matrices
 * - Proper response structure and data validation
 * - Error handling and edge cases
 */

const request = require('supertest');
const { createApp } = require('../../src/server/api/server');
const { createTestContainer } = require('../../src/core/services');
const {
  createTestSetup,
  ERROR_MESSAGES,
  STATUS_CODES
} = require('../../agents/structuralize-api-contract');

/**
 * Test data for complex chemical matrices
 * Each represents a real-world substance with multiple chemical compounds
 */
const COMPLEX_CHEMICAL_MATRICES = {
  cacao: {
    input: 'cacao',
    chemicals: [
      { name: "Theobromine", smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C.O" },
      { name: "Caffeine", smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" },
      { name: "Epicatechin", smiles: "C1C(C(OC2=CC(=CC(=C21)O)O)C3=CC(=C(C(=C3)O)O)O)O" },
      { name: "Catechin", smiles: "C1C(C(OC2=CC(=CC(=C21)O)O)C3=CC(=C(C(=C3)O)O)O)O" },
      { name: "Anandamide", smiles: "CCCCCCCCC(=O)NC(CO)C" }
    ],
    reason: "Cacao contains methylxanthines (theobromine, caffeine), flavonoids (epicatechin, catechin), and bioactive lipids (anandamide)"
  },

  coffee: {
    input: 'coffee',
    chemicals: [
      { name: "Caffeine", smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" },
      { name: "Chlorogenic Acid", smiles: "C1C(C(C(CC1(C(=O)O)O)OC(=O)C=CC2=CC(=C(C=C2)O)O)O)O" },
      { name: "Caffeic Acid", smiles: "C1=CC(=C(C=C1C=CC(=O)O)O)O" },
      { name: "Trigonelline", smiles: "C[N+]1=CC=CC(=C1)C(=O)O" }
    ],
    reason: "Coffee contains alkaloids (caffeine, trigonelline) and phenolic acids (chlorogenic, caffeic)"
  },

  wine: {
    input: 'wine',
    chemicals: [
      { name: "Ethanol", smiles: "CCO" },
      { name: "Resveratrol", smiles: "C1=CC(=CC=C1C=CC2=CC(=CC(=C2)O)O)O" },
      { name: "Quercetin", smiles: "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O" },
      { name: "Tartaric Acid", smiles: "C(C(C(=O)O)O)(C(=O)O)O" }
    ],
    reason: "Wine contains ethanol, polyphenols (resveratrol, quercetin), and organic acids (tartaric acid)"
  },

  // Example: Adding a new complex matrix is trivial - just add to this object!
  tea: {
    input: 'tea',
    chemicals: [
      { name: "Catechin", smiles: "C1C(C(OC2=CC(=CC(=C21)O)O)C3=CC(=C(C(=C3)O)O)O)O" },
      { name: "Epigallocatechin Gallate", smiles: "C1C(C(OC2=CC(=CC(=C21)O)O)C3=CC(=C(C(=C3)O)O)O)OC(=O)C4=CC(=C(C(=C4)O)O)O" },
      { name: "Caffeine", smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" },
      { name: "Theanine", smiles: "C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N" }
    ],
    reason: "Tea contains catechins (catechin, EGCG), caffeine, and amino acids (theanine)"
  }
};

describe('Complex Chemical Analysis Integration Test', () => {
  let app;
  let container;
  let mocks;
  const testSetup = createTestSetup(createTestContainer, createApp);

  beforeEach(async () => {
    const setup = await testSetup.setup();
    app = setup.app;
    container = setup.container;
    mocks = setup.mocks;
  });

  afterEach(() => {
    testSetup.teardown({ container });
  });

  describe('Complex Chemical Matrix Analysis', () => {
    // Parameterized test for all complex chemical matrices
    Object.entries(COMPLEX_CHEMICAL_MATRICES).forEach(([matrixName, matrixData]) => {
      test(`should successfully analyze ${matrixName} and return multiple chemical compounds`, async () => {
        // Mock AI response for the specific chemical matrix
        mocks.openaiClient.chat.completions.create.mockResolvedValue({
          choices: [{
            message: {
              content: JSON.stringify({
                object: matrixData.input,
                chemicals: matrixData.chemicals,
                reason: matrixData.reason
              })
            }
          }]
        });

        const response = await request(app)
          .post('/api/structuralize')
          .send({
            text: matrixData.input,
            lookupMode: 'GPT-5'
          })
          .expect(STATUS_CODES.SUCCESS);

        expect(response.body).toHaveProperty('object', matrixData.input);
        expect(response.body).toHaveProperty('chemicals');
        expect(Array.isArray(response.body.chemicals)).toBe(true);
        expect(response.body.chemicals.length).toBeGreaterThan(1); // Multiple compounds

        // Validate each chemical has required properties
        response.body.chemicals.forEach((chemical, index) => {
          expect(chemical).toHaveProperty('name');
          expect(chemical).toHaveProperty('smiles');
          expect(typeof chemical.name).toBe('string');
          expect(typeof chemical.smiles).toBe('string');
          expect(chemical.smiles.length).toBeGreaterThan(0);

          // Verify the specific chemicals expected for this matrix
          expect(chemical.name).toBe(matrixData.chemicals[index].name);
          expect(chemical.smiles).toBe(matrixData.chemicals[index].smiles);
        });

        // Verify AI service was called correctly with the right input
        expect(mocks.openaiClient.chat.completions.create).toHaveBeenCalledWith(
          expect.objectContaining({
            messages: expect.arrayContaining([
              expect.objectContaining({
                role: 'user',
                content: expect.stringContaining(matrixData.input)
              })
            ])
          })
        );
      });
    });

    test('should handle input variations and synonyms', async () => {
      const inputVariations = [
        { input: 'cacao', synonym: 'cocoa' },
        { input: 'coffee', synonym: 'java' },
        { input: 'wine', synonym: 'grape wine' }
      ];

      for (const { input, synonym } of inputVariations) {
        const matrixData = COMPLEX_CHEMICAL_MATRICES[input];

        mocks.openaiClient.chat.completions.create.mockResolvedValue({
          choices: [{
            message: {
              content: JSON.stringify({
                object: synonym,
                chemicals: matrixData.chemicals.slice(0, 2), // Just first 2 compounds for variation test
                reason: `${synonym} contains ${matrixData.chemicals[0].name} and ${matrixData.chemicals[1].name}`
              })
            }
          }]
        });

        const response = await request(app)
          .post('/api/structuralize')
          .send({ text: synonym, lookupMode: 'GPT-5' })
          .expect(STATUS_CODES.SUCCESS);

        expect(response.body.object).toBe(synonym);
        expect(response.body.chemicals.length).toBeGreaterThan(0);
      }
    });

    test('should validate SMILES notation format for all complex matrices', async () => {
      // Test SMILES validation for each chemical matrix
      for (const [matrixName, matrixData] of Object.entries(COMPLEX_CHEMICAL_MATRICES)) {
        mocks.openaiClient.chat.completions.create.mockResolvedValue({
          choices: [{
            message: {
              content: JSON.stringify({
                object: matrixData.input,
                chemicals: matrixData.chemicals,
                reason: matrixData.reason
              })
            }
          }]
        });

        const response = await request(app)
          .post('/api/structuralize')
          .send({ text: matrixData.input, lookupMode: 'GPT-5' })
          .expect(STATUS_CODES.SUCCESS);

        response.body.chemicals.forEach(chemical => {
          // Basic SMILES validation - should contain common atoms and no invalid characters
          expect(chemical.smiles).toMatch(/^[A-Za-z0-9\[\]\(\)=\-#+\.]+$/);
          expect(chemical.smiles.length).toBeGreaterThan(2);

          // Verify SMILES matches expected values from test data
          const expectedChemical = matrixData.chemicals.find(c => c.name === chemical.name);
          expect(expectedChemical).toBeDefined();
          expect(chemical.smiles).toBe(expectedChemical.smiles);
        });
      }
    });

    test('should handle AI service timeout gracefully for complex matrices', async () => {
      // Test timeout handling for each matrix type
      for (const [matrixName, matrixData] of Object.entries(COMPLEX_CHEMICAL_MATRICES)) {
        // Mock AI service timeout for each matrix
        mocks.openaiClient.chat.completions.create.mockRejectedValue(
          new Error(`OpenAI API timeout for ${matrixData.input}`)
        );

        const response = await request(app)
          .post('/api/structuralize')
          .send({ text: matrixData.input, lookupMode: 'GPT-5' })
          .expect(STATUS_CODES.INTERNAL_ERROR);

        expect(response.body).toHaveProperty('error');
        expect(response.body.error).toContain('failed');
        expect(response.body).toHaveProperty('code');
      }
    });

    test('should reject empty cacao input', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: '', lookupMode: 'GPT-5' })
        .expect(STATUS_CODES.BAD_REQUEST);

      expect(response.body.error).toContain(ERROR_MESSAGES.MISSING_TEXT);
    });

    test('should reject missing text parameter', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ lookupMode: 'GPT-5' })
        .expect(STATUS_CODES.BAD_REQUEST);

      expect(response.body.error).toContain(ERROR_MESSAGES.MISSING_TEXT);
    });
  });

  describe('Cacao Analysis Performance', () => {
    test('should complete analysis within reasonable time for all matrices', async () => {
      // Test performance for each complex matrix
      for (const [matrixName, matrixData] of Object.entries(COMPLEX_CHEMICAL_MATRICES)) {
        mocks.openaiClient.chat.completions.create.mockResolvedValue({
          choices: [{
            message: {
              content: JSON.stringify({
                object: matrixData.input,
                chemicals: matrixData.chemicals.slice(0, 2), // Just first 2 for performance test
                reason: `${matrixData.input} analysis`
              })
            }
          }]
        });

        const startTime = Date.now();

        await request(app)
          .post('/api/structuralize')
          .send({ text: matrixData.input, lookupMode: 'GPT-5' })
          .expect(STATUS_CODES.SUCCESS);

        const duration = Date.now() - startTime;
        expect(duration).toBeLessThan(5000); // Should complete within 5 seconds

        console.log(`✅ ${matrixName} analysis completed in ${duration}ms`);
      }
    });
  });
});