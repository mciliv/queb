/**
 * Server Integration Tests with Dependency Injection
 * 
 * These tests demonstrate how DI makes integration testing easier by:
 * 1. Allowing partial mocking (mock only external services)
 * 2. Testing with different configurations easily
 * 3. Creating isolated test environments
 * 4. Avoiding side effects between tests
 */

const request = require('supertest');
const { createApp } = require('../../src/server/api/server');
const { createTestContainer } = require('../../src/core/services');

describe('Server Integration Tests with DI', () => {
  describe('API Endpoints', () => {
    let app;
    let container;
    let mocks;

    beforeEach(async () => {
      // Create test container with partial mocks
      // We mock only external services, leaving internal logic real
      mocks = {
        openaiClient: {
          chat: {
            completions: {
              create: jest.fn()
            }
          }
        },
        // Real services will be used for:
        // - molecularProcessor
        // - errorHandler
        // - promptEngine
        // This gives us confidence that integration works correctly
      };

      container = createTestContainer(mocks);
      app = await createApp(container);
    });

    afterEach(() => {
      // Clean up is automatic - no global state to worry about!
      container.clear();
    });

    describe('POST /api/structuralize', () => {
      it('should analyze text and return chemicals', async () => {
        // Setup AI mock
        mocks.openaiClient.chat.completions.create.mockResolvedValue({
          choices: [{
            message: {
              content: JSON.stringify({
                object: 'water',
                chemicals: [
                  { name: 'Water', smiles: 'O' }
                ]
              })
            }
          }]
        });

        // Make request
        const response = await request(app)
          .post('/api/structuralize')
          .send({ text: 'water' })
          .expect(200);

        // Verify response
        expect(response.body).toMatchObject({
          object: 'water',
          chemicals: expect.arrayContaining([
            expect.objectContaining({
              name: 'Water',
              smiles: 'O'
            })
          ])
        });
      });

      it('should validate input', async () => {
        const response = await request(app)
          .post('/api/structuralize')
          .send({ text: 123 }) // Invalid type
          .expect(400);

        expect(response.body.error).toContain('invalid');
      });

      it('should handle AI service errors gracefully', async () => {
        // Setup AI to fail
        mocks.openaiClient.chat.completions.create.mockRejectedValue(
          new Error('OpenAI API error')
        );

        const response = await request(app)
          .post('/api/structuralize')
          .send({ text: 'test' })
          .expect(500);

        expect(response.body.error).toBeDefined();
      });
    });

    describe('POST /api/structuralize-image', () => {
      it('should analyze image with coordinates', async () => {
        // Mock detection and prediction
        mocks.openaiClient.chat.completions.create
          .mockResolvedValueOnce({
            choices: [{
              message: {
                content: JSON.stringify({
                  object: 'coffee cup',
                  recommendedBox: { x: 10, y: 10, width: 100, height: 100 }
                })
              }
            }]
          })
          .mockResolvedValueOnce({
            choices: [{
              message: {
                content: JSON.stringify({
                  object: 'coffee cup',
                  chemicals: [{ name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' }]
                })
              }
            }]
          });

        const response = await request(app)
          .post('/api/structuralize-image')
          .send({
            imageBase64: 'base64-test-data',
            x: 50,
            y: 50
          })
          .expect(200);

        expect(response.body.object).toBe('coffee cup');
        expect(response.body.recommendedBox).toBeDefined();
      });
    });
  });

  describe('Configuration Variations', () => {
    it('should work without AI when configured', async () => {
      // Create container with no AI
      const container = createTestContainer({
        openaiClient: null,
        config: {
          get: jest.fn((key) => {
            if (key === 'openai.apiKey') return null;
            if (key === 'database.enabled') return false;
            return null;
          }),
          isProduction: jest.fn(() => false),
          isDevelopment: jest.fn(() => true),
          validate: jest.fn()
        }
      });

      const app = await createApp(container);

      // Text prediction should fail without AI
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: 'test' })
        .expect(500);

      expect(response.body.error).toContain('AI');
    });

    it('should enable user routes when database is configured', async () => {
      // Create container with database
      const mockUserService = {
        getUserByDeviceToken: jest.fn(),
        createUser: jest.fn(),
        updateUser: jest.fn()
      };

      const container = createTestContainer({
        config: {
          get: jest.fn((key) => {
            if (key === 'database.enabled') return true;
            return null;
          }),
          isProduction: jest.fn(() => false),
          isDevelopment: jest.fn(() => true),
          validate: jest.fn()
        },
        userService: mockUserService
      });

      const app = await createApp(container);

      // User endpoint should exist
      mockUserService.getUserByDeviceToken.mockResolvedValue({
        id: 1,
        deviceToken: 'test-token'
      });

      const response = await request(app)
        .get('/api/user/test-token')
        .expect(200);

      expect(response.body.deviceToken).toBe('test-token');
    });
  });

  describe('Performance and Caching', () => {
    it('should use cache when available', async () => {
      // Create container with cache
      const cacheHits = new Map();
      const mockCache = {
        get: jest.fn(key => cacheHits.get(key)),
        set: jest.fn((key, value) => cacheHits.set(key, value))
      };

      const container = createTestContainer({
        cache: mockCache,
        openaiClient: {
          chat: {
            completions: {
              create: jest.fn().mockResolvedValue({
                choices: [{
                  message: {
                    content: JSON.stringify({
                      object: 'water',
                      chemicals: [{ name: 'Water', smiles: 'O' }]
                    })
                  }
                }]
              })
            }
          }
        }
      });

      const app = await createApp(container);

      // First request - cache miss
      await request(app)
        .post('/api/structuralize')
        .send({ text: 'water' })
        .expect(200);

      expect(mockCache.set).toHaveBeenCalled();

      // Second request - should hit cache
      await request(app)
        .post('/api/structuralize')
        .send({ text: 'water' })
        .expect(200);

      // AI should only be called once due to caching
      const aiClient = await container.get('openaiClient');
      expect(aiClient.chat.completions.create).toHaveBeenCalledTimes(1);
    });
  });

  describe('Error Scenarios', () => {
    it('should handle database connection failures', async () => {
      const container = createTestContainer({
        database: {
          initialize: jest.fn().mockRejectedValue(new Error('Connection failed'))
        },
        config: {
          get: jest.fn((key) => {
            if (key === 'database.enabled') return true;
            return null;
          }),
          validate: jest.fn()
        }
      });

      // Server should still start even if DB fails
      const app = await createApp(container);
      
      const response = await request(app)
        .get('/api/health')
        .expect(200);

      expect(response.body.status).toBe('ok');
    });
  });
});

/**
 * Benefits demonstrated:
 * 
 * 1. **Isolated Tests**: Each test gets its own container and app instance
 * 2. **Flexible Mocking**: Can mock just what we need, keep rest real
 * 3. **Configuration Testing**: Easy to test different configurations
 * 4. **No Global State**: No need to reset modules or worry about test order
 * 5. **Clear Dependencies**: Tests clearly show what dependencies are used
 * 6. **Fast Tests**: Can mock slow external services while testing real logic
 * 
 * Without DI, we'd need:
 * - Complex jest.mock() setups
 * - Module cache manipulation
 * - Careful test ordering
 * - Environment variable manipulation
 * - Singleton reset logic
 */
