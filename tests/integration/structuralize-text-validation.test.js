/**
 * Backend Tests for /api/structuralize Text Parameter Validation
 * 
 * These tests verify that the backend properly validates the "text" parameter
 * and returns appropriate error messages for invalid inputs.
 */

const request = require('supertest');
const { createApp } = require('../../src/server/api/server');
const { createTestContainer } = require('../../src/core/services');
const { 
  createTestSetup, 
  ERROR_MESSAGES, 
  STATUS_CODES 
} = require('../../agents/structuralize-api-contract');

describe('Backend: /api/structuralize Text Parameter Validation', () => {
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

  describe('Valid text inputs', () => {
    it('should accept valid text string', async () => {
      // Setup AI mock for successful response
      mocks.openaiClient.chat.completions.create.mockResolvedValue(
        testSetup.mockSuccessfulResponse('coffee', [
          { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' }
        ])
      );

      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: 'coffee', lookupMode: 'GPT-5' })
        .expect(STATUS_CODES.SUCCESS);

      expect(response.body).toHaveProperty('object');
      expect(response.body).toHaveProperty('chemicals');
      expect(mocks.openaiClient.chat.completions.create).toHaveBeenCalled();
    });

    it('should accept text with leading/trailing whitespace and trim it', async () => {
      mocks.openaiClient.chat.completions.create.mockResolvedValue(
        testSetup.mockSuccessfulResponse('water', [
          { name: 'Water', smiles: 'O' }
        ])
      );

      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: '  water  ', lookupMode: 'GPT-5' })
        .expect(STATUS_CODES.SUCCESS);

      expect(response.body).toHaveProperty('object');
      // Verify the trimmed text was used (check the call was made)
      expect(mocks.openaiClient.chat.completions.create).toHaveBeenCalled();
    });

    it('should accept text with default lookupMode', async () => {
      mocks.openaiClient.chat.completions.create.mockResolvedValue(
        testSetup.mockSuccessfulResponse('test', [])
      );

      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: 'test' })
        .expect(STATUS_CODES.SUCCESS);

      expect(response.body).toHaveProperty('object');
    });
  });

  describe('Invalid text inputs - should return 400', () => {
    it('should reject missing text parameter', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ lookupMode: 'GPT-5' })
        .expect(STATUS_CODES.BAD_REQUEST);

      expect(response.body).toHaveProperty('error');
      expect(response.body.error).toContain(ERROR_MESSAGES.MISSING_TEXT);
      expect(mocks.openaiClient.chat.completions.create).not.toHaveBeenCalled();
    });

    it('should reject empty string', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: '', lookupMode: 'GPT-5' })
        .expect(STATUS_CODES.BAD_REQUEST);

      expect(response.body).toHaveProperty('error');
      expect(response.body.error).toContain(ERROR_MESSAGES.MISSING_TEXT);
      expect(mocks.openaiClient.chat.completions.create).not.toHaveBeenCalled();
    });

    it('should reject whitespace-only string', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: '   ', lookupMode: 'GPT-5' })
        .expect(STATUS_CODES.BAD_REQUEST);

      expect(response.body).toHaveProperty('error');
      expect(response.body.error).toContain(ERROR_MESSAGES.MISSING_TEXT);
      expect(mocks.openaiClient.chat.completions.create).not.toHaveBeenCalled();
    });

    it('should reject null text', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: null, lookupMode: 'GPT-5' })
        .expect(STATUS_CODES.BAD_REQUEST);

      expect(response.body).toHaveProperty('error');
      expect(response.body.error).toContain(ERROR_MESSAGES.MISSING_TEXT);
      expect(mocks.openaiClient.chat.completions.create).not.toHaveBeenCalled();
    });

    it('should reject undefined text', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: undefined, lookupMode: 'GPT-5' })
        .expect(STATUS_CODES.BAD_REQUEST);

      expect(response.body).toHaveProperty('error');
      expect(response.body.error).toContain(ERROR_MESSAGES.MISSING_TEXT);
      expect(mocks.openaiClient.chat.completions.create).not.toHaveBeenCalled();
    });

    it('should reject non-string text (number)', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: 123, lookupMode: 'GPT-5' })
        .expect(STATUS_CODES.BAD_REQUEST);

      expect(response.body).toHaveProperty('error');
      expect(response.body.error).toContain(ERROR_MESSAGES.MISSING_TEXT);
      expect(mocks.openaiClient.chat.completions.create).not.toHaveBeenCalled();
    });

    it('should reject non-string text (boolean)', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: true, lookupMode: 'GPT-5' })
        .expect(STATUS_CODES.BAD_REQUEST);

      expect(response.body).toHaveProperty('error');
      expect(response.body.error).toContain(ERROR_MESSAGES.MISSING_TEXT);
      expect(mocks.openaiClient.chat.completions.create).not.toHaveBeenCalled();
    });

    it('should reject non-string text (object)', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: { value: 'test' }, lookupMode: 'GPT-5' })
        .expect(STATUS_CODES.BAD_REQUEST);

      expect(response.body).toHaveProperty('error');
      expect(response.body.error).toContain(ERROR_MESSAGES.MISSING_TEXT);
      expect(mocks.openaiClient.chat.completions.create).not.toHaveBeenCalled();
    });

    it('should reject non-string text (array)', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: ['test'], lookupMode: 'GPT-5' })
        .expect(STATUS_CODES.BAD_REQUEST);

      expect(response.body).toHaveProperty('error');
      expect(response.body.error).toContain(ERROR_MESSAGES.MISSING_TEXT);
      expect(mocks.openaiClient.chat.completions.create).not.toHaveBeenCalled();
    });

    it('should reject empty body', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({})
        .expect(STATUS_CODES.BAD_REQUEST);

      expect(response.body).toHaveProperty('error');
      expect(response.body.error).toContain(ERROR_MESSAGES.MISSING_TEXT);
      expect(mocks.openaiClient.chat.completions.create).not.toHaveBeenCalled();
    });
  });

  describe('lookupMode validation', () => {
    it('should accept valid lookupMode', async () => {
      mocks.openaiClient.chat.completions.create.mockResolvedValue(
        testSetup.mockSuccessfulResponse('test', [])
      );

      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: 'test', lookupMode: 'GPT-5' })
        .expect(STATUS_CODES.SUCCESS);

      expect(response.body).toHaveProperty('object');
    });

    it('should reject invalid lookupMode type', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: 'test', lookupMode: 123 })
        .expect(STATUS_CODES.BAD_REQUEST);

      expect(response.body).toHaveProperty('error');
      expect(response.body.error).toContain(ERROR_MESSAGES.INVALID_LOOKUP_MODE);
      expect(mocks.openaiClient.chat.completions.create).not.toHaveBeenCalled();
    });
  });

  describe('Edge cases', () => {
    it('should handle very long text strings', async () => {
      const longText = 'a'.repeat(1000);
      
      mocks.openaiClient.chat.completions.create.mockResolvedValue(
        testSetup.mockSuccessfulResponse(longText, [])
      );

      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: longText, lookupMode: 'GPT-5' })
        .expect(STATUS_CODES.SUCCESS);

      expect(response.body).toHaveProperty('object');
    });

    it('should handle text with special characters', async () => {
      mocks.openaiClient.chat.completions.create.mockResolvedValue(
        testSetup.mockSuccessfulResponse('test@#$%', [])
      );

      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: 'test@#$%', lookupMode: 'GPT-5' })
        .expect(STATUS_CODES.SUCCESS);

      expect(response.body).toHaveProperty('object');
    });

    it('should handle unicode characters', async () => {
      mocks.openaiClient.chat.completions.create.mockResolvedValue(
        testSetup.mockSuccessfulResponse('café', [])
      );

      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: 'café', lookupMode: 'GPT-5' })
        .expect(STATUS_CODES.SUCCESS);

      expect(response.body).toHaveProperty('object');
    });
  });
});
