const request = require('supertest');
const { createApp } = require('../../../src/server/api/app');
const { createContainer } = require('../../../src/core/services');

describe('POST /api/log-error - Clean Architecture Integration', () => {
  let app;
  let container;
  let mockLogErrorPort;

  beforeAll(async () => {
    // SETUP CONTAINER: Wire up the clean architecture
    container = createContainer();

    // MOCK PORT: Create a mock that implements the port interface
    mockLogErrorPort = {
      log: jest.fn()
    };

    // Override the use-case directly to use our mock port
    container.register('logErrorUseCase', () => {
      const LogErrorUseCase = require('../../../src/server/services/log-error-use-case');
      return new LogErrorUseCase({ logErrorPort: mockLogErrorPort });
    });

    // CREATE APP: With our mocked container
    app = await createApp({
      config: { get: () => ({}) }, // Minimal config
      logger: { info: jest.fn(), error: jest.fn() }, // Mock logger
      container
    });
  });

  beforeEach(() => {
    jest.clearAllMocks();
  });

  test('should return 200 and delegate to use-case', async () => {
    // INTEGRATION TEST: Full HTTP → domain flow
    const response = await request(app)
      .post('/api/log-error')
      .send({
        type: 'error',
        message: 'Test error from frontend',
        source: 'frontend',
        location: 'App.jsx:123'
      })
      .expect(200);

    expect(response.body).toEqual({ received: true });

    // VERIFY PORT CALLED: Domain logic executed correctly
    expect(mockLogErrorPort.log).toHaveBeenCalledWith({
      level: 'error',
      message: '[frontend] error: Test error from frontend',
      metadata: expect.objectContaining({
        type: 'error',
        message: 'Test error from frontend',
        source: 'frontend',
        location: 'App.jsx:123'
      })
    });
  });

  test('should handle use-case errors gracefully', async () => {
    // ERROR TEST: Verify HTTP error handling when domain fails
    mockLogErrorPort.log.mockRejectedValue(new Error('Logger failed'));

    const response = await request(app)
      .post('/api/log-error')
      .send({
        type: 'error',
        message: 'Will fail'
      })
      .expect(500);

    expect(response.body).toEqual({ error: 'Failed to log error' });
  });
});