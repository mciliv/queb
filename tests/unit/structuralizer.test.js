
const Structuralizer = require('../../src/server/services/Structuralizer');

describe('Structuralizer with Dependency Injection', () => {
  // Mock factories for easy test setup
  const createMocks = () => ({
    aiService: {
      callAPI: jest.fn()
    },
    molecularProcessor: {
      generateSDF: jest.fn(),
      processSmiles: jest.fn()
    },
    nameResolver: {
      resolveName: jest.fn()
    },
    promptEngine: {
      generateChemicalPrompt: jest.fn(),
      generateDetectionPrompt: jest.fn(),
      validateResponse: jest.fn(),
      repairJSON: jest.fn()
    },
    errorHandler: {
      handle: jest.fn(error => ({
        message: error.message,
        code: 'TEST_ERROR',
        statusCode: 500
      }))
    },
    logger: {
      info: jest.fn(),
      warn: jest.fn(),
      error: jest.fn()
    },
    cache: {
      get: jest.fn(),
      set: jest.fn()
    },
    config: {
      cacheEnabled: true, // Enable caching for tests
      get: jest.fn((key) => {
        const configMap = {
          'openai.model': 'gpt-3.5-turbo',
          'openai.timeout': 30000,
          'openai.useCompletionAPI': false, // Force chat completions for tests
          'cache.maxSize': 100,
          'cache.ttl': 300000
        };
        return configMap[key];
      })
    },
    foodDbService: {
      getCompounds: jest.fn()
    },
    wikidataResolver: {
      findChemicals: jest.fn()
    }
  });

  describe('Construction and Validation', () => {
    it('should validate required dependencies', () => {
      expect(() => new Structuralizer({})).toThrow('Missing required dependencies');
    });

    it('should accept all dependencies', () => {
      const mocks = createMocks();
      const structuralizer = new Structuralizer(mocks);
      
      expect(structuralizer.aiService).toBe(mocks.aiService);
      expect(structuralizer.logger).toBe(mocks.logger);
    });
  });

  describe('Text Prediction', () => {
    let structuralizer;
    let mocks;

    beforeEach(() => {
      mocks = createMocks();
      structuralizer = new Structuralizer(mocks);
    });

    it('should analyze text using AI', async () => {
      // Setup mocks
      mocks.promptEngine.generateChemicalPrompt.mockReturnValue('Analyze coffee');
      mocks.promptEngine.validateResponse.mockReturnValue(true);

      // Mock OpenAI service - returns parsed JSON from AI
      mocks.aiService.callAPI.mockResolvedValue({
        object: 'coffee',
        chemicals: [
          { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' }
        ]
      });
      mocks.molecularProcessor.generateSDF.mockResolvedValue('/sdf/caffeine.sdf');

      // Execute
      const result = await structuralizer.chemicals({
        object: 'coffee',
        lookupMode: 'ai'
      });

      // Verify
      expect(result.object).toBe('coffee');
      expect(result.chemicals).toHaveLength(1);
      // Structuralizer returns the AI result; SDF generation happens via /api/generate-sdfs
      expect(result.chemicals[0]).toEqual({
        name: 'Caffeine',
        smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
      });

      // Verify calls
      expect(mocks.promptEngine.generateChemicalPrompt).toHaveBeenCalledWith(
        'coffee',
        { includeReason: true }
      );
      expect(mocks.molecularProcessor.generateSDF).not.toHaveBeenCalled();
    });

    it('should not fail validation when some chemicals have missing SMILES (resolve by name later)', async () => {
      // Use real validation rules for this regression test
      const realPromptEngine = require('../../src/core/PromptEngine');
      mocks.promptEngine.generateChemicalPrompt.mockReturnValue('Analyze cacao');
      mocks.promptEngine.validateResponse = realPromptEngine.validateResponse.bind(realPromptEngine);

      // AI returns some compounds without smiles (this used to cause "Load failed")
      mocks.aiService.callAPI.mockResolvedValue({
        object: 'cacao',
        chemicals: [
          { name: 'Theobromine', smiles: 'CN1C=NC2=C1C(=O)NC(=O)N2' },
          { name: 'Caffeine', smiles: null },
        ]
      });

      const result = await structuralizer.chemicals({
        object: 'cacao',
        lookupMode: 'ai'
      });

      expect(result.object).toBe('cacao');
      expect(result.chemicals).toHaveLength(2);
      expect(result.chemicals[0].name).toBe('Theobromine');
      expect(result.chemicals[1].name).toBe('Caffeine');
      expect(mocks.molecularProcessor.generateSDF).not.toHaveBeenCalled();
    });

    it('should use cache when available', async () => {
      // Setup cache hit
      const cachedResult = {
        object: 'water',
        chemicals: [{ name: 'Water', smiles: 'O', sdfPath: '/sdf/water.sdf' }]
      };
      mocks.cache.get.mockResolvedValue(cachedResult);

      // Execute
      const result = await structuralizer.chemicals({
        object: 'water',
        lookupMode: 'ai'
      });

      // Verify
      expect(result).toEqual(cachedResult);
      expect(mocks.aiService.callAPI).not.toHaveBeenCalled();
      expect(mocks.cache.get).toHaveBeenCalled();
    });

    it.skip('should handle FoodDB lookups', async () => {
      // TODO: Implement FoodDB lookup mode
      // Setup mocks
      mocks.foodDbService.getCompounds.mockResolvedValue([
        { name: 'Glucose', moldb_smiles: 'C(C1C(C(C(C(O1)O)O)O)O)O' },
        { name: 'Fructose', moldb_smiles: 'C(C1C(C(C(O1)(CO)O)O)O)O' }
      ]);
      mocks.molecularProcessor.generateSDF
        .mockResolvedValueOnce('/sdf/glucose.sdf')
        .mockResolvedValueOnce('/sdf/fructose.sdf');

      // Execute
      const result = await structuralizer.chemicals({
        object: 'apple',
        lookupMode: 'foodb'
      });

      // Verify
      expect(result.chemicals).toHaveLength(2);
      expect(mocks.foodDbService.getCompounds).toHaveBeenCalledWith('apple');
      expect(mocks.openAIService.callAPI).not.toHaveBeenCalled();
    });

    it.skip('should fallback to AI when FoodDB returns no results', async () => {
      // TODO: Implement FoodDB lookup mode
      // Setup mocks
      mocks.foodDbService.getCompounds.mockResolvedValue([]);
      mocks.promptEngine.generateChemicalPrompt.mockReturnValue('Analyze rare compound');
      mocks.promptEngine.validateResponse.mockReturnValue(true);
      mocks.openAIService.callAPI.mockResolvedValue({
        choices: [{
          message: {
            content: JSON.stringify({
              object: 'rare compound',
              chemicals: [{ name: 'Unknown', smiles: null }]
            })
          }
        }]
      });

      mocks.aiClient.completions.create.mockResolvedValue({
        choices: [{
          text: JSON.stringify({
            object: 'rare compound',
            chemicals: [{ name: 'Unknown', smiles: null }]
          })
        }]
      });

      // Execute
      const result = await structuralizer.chemicals({
        object: 'rare compound',
        lookupMode: 'foodb'
      });

      // Verify fallback
      expect(mocks.foodDbService.getCompounds).toHaveBeenCalled();
      expect(mocks.openAIService.callAPI).toHaveBeenCalled();
    });
  });

  describe.skip('Image Prediction', () => {
    let structuralizer;
    let mocks;

    beforeEach(() => {
      mocks = createMocks();
      structuralizer = new Structuralizer(mocks);
    });

    it('should detect object in image and analyze', async () => {
      // Setup mocks
      mocks.promptEngine.generateDetectionPrompt.mockReturnValue('Detect object at (100, 200)');
      mocks.promptEngine.generateChemicalPrompt.mockReturnValue('Analyze detected object');
      mocks.promptEngine.validateResponse.mockReturnValue(true);
      
      // Mock OpenAI service calls
      mocks.aiService.callAPI
        .mockResolvedValueOnce({
          object: 'coffee cup',
          recommendedBox: { x: 50, y: 150, width: 100, height: 100 }
        })
        .mockResolvedValueOnce({
          object: 'coffee cup',
          chemicals: [{ name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' }]
        });

      mocks.molecularProcessor.generateSDF.mockResolvedValue('/sdf/caffeine.sdf');

      // Execute
      const result = await structuralizer.chemicals({
        imageBase64: 'base64-image-data',
        x: 100,
        y: 200
      });

      // Verify
      expect(result.object).toBe('coffee cup');
      expect(result.recommendedBox).toEqual({ x: 50, y: 150, width: 100, height: 100 });
      expect(result.chemicals).toHaveLength(1);
      
      // Verify two AI calls (detection + analysis)
      expect(mocks.aiService.callAPI).toHaveBeenCalledTimes(2);
    });
  });

  describe('Error Handling', () => {
    let structuralizer;
    let mocks;

    beforeEach(() => {
      mocks = createMocks();
      structuralizer = new Structuralizer(mocks);
    });

    it.skip('should handle AI timeout errors with retry', async () => {
      // TODO: Implement retry logic in Structuralizer
      // Setup mocks
      mocks.promptEngine.generateChemicalPrompt.mockReturnValue('Analyze');

      // Mock OpenAI service with retry behavior
      mocks.openAIService.callAPI
        .mockRejectedValueOnce(new Error('timeout'))
        .mockResolvedValueOnce({
          object: 'water',
          chemicals: [{ name: 'Water', smiles: 'O' }]
        });

      mocks.promptEngine.validateResponse.mockReturnValue(true);
      mocks.molecularProcessor.generateSDF.mockResolvedValue('/sdf/water.sdf');

      // Execute
      const result = await structuralizer.chemicals({
        object: 'water',
        lookupMode: 'ai'
      });

      // Verify retry worked
      expect(result.object).toBe('water');
      expect(mocks.openAIService.callAPI).toHaveBeenCalledTimes(2);
      expect(mocks.logger.warn).toHaveBeenCalledWith(
        'AI call failed (attempt 1)',
        { error: 'timeout' }
      );
    });

    it('should handle invalid AI responses', async () => {
      // Setup mocks
      mocks.promptEngine.generateChemicalPrompt.mockReturnValue('Analyze');
      // Mock the OpenAI service to throw an error for invalid responses
      mocks.aiService.callAPI.mockRejectedValue(
        new Error('AI response parse failed: Unexpected token \'I\'')
      );

      // Execute and expect error
      await expect(structuralizer.chemicals({
        object: 'test',
        lookupMode: 'ai'
      })).rejects.toThrow('AI response parse failed');
    });
  });

  describe.skip('Structure Generation', () => {
    let structuralizer;
    let mocks;

    beforeEach(() => {
      mocks = createMocks();
      structuralizer = new Structuralizer(mocks);
    });

    it('should handle molecules without SMILES by name resolution', async () => {
      // Setup mocks
      mocks.promptEngine.generateChemicalPrompt.mockReturnValue('Analyze');
      mocks.promptEngine.validateResponse.mockReturnValue(true);
      mocks.aiService.callAPI.mockResolvedValue({
        object: 'aspirin',
        chemicals: [{ name: 'Aspirin', smiles: null }] // No SMILES
      });

      // Mock name resolution
      mocks.nameResolver.resolveName.mockResolvedValue({
        name: 'Aspirin',
        smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O',
        iupac: 'acetylsalicylic acid'
      });

      mocks.molecularProcessor.generateSDF.mockResolvedValue('/sdf/aspirin.sdf');

      // Execute
      const result = await structuralizer.chemicals({
        object: 'aspirin',
        lookupMode: 'ai'
      });

      // Verify
      expect(result.chemicals[0].status).toBe('ok');
      expect(result.chemicals[0].sdfPath).toBe('/sdf/aspirin.sdf');
      expect(mocks.nameResolver.resolveName).toHaveBeenCalledWith('Aspirin');
    });

    it('should handle structure generation failures gracefully', async () => {
      // Setup mocks
      mocks.promptEngine.generateChemicalPrompt.mockReturnValue('Analyze');
      mocks.promptEngine.validateResponse.mockReturnValue(true);
      mocks.aiService.callAPI.mockResolvedValue({
        object: 'complex',
        chemicals: [
          { name: 'Valid', smiles: 'CCO' },
          { name: 'Invalid', smiles: 'INVALID_SMILES' }
        ]
      });

      mocks.molecularProcessor.generateSDF
        .mockResolvedValueOnce('/sdf/valid.sdf')
        .mockRejectedValueOnce(new Error('Invalid SMILES'));

      // Execute
      const result = await structuralizer.chemicals({
        object: 'complex',
        lookupMode: 'ai'
      });

      // Verify partial success
      expect(result.chemicals).toHaveLength(2);
      expect(result.chemicals[0].status).toBe('ok');
      expect(result.chemicals[1].status).toBe('error');
      expect(mocks.logger.warn).toHaveBeenCalled();
    });
  });
});
