/**
 * Unit tests for text structuralize endpoint validation
 * Tests input validation rules without full server setup
 */

const request = require('supertest');
const express = require('express');

describe('POST /api/structuralize - Input Validation', () => {
  let app;
  let mockStructuralizer;

  beforeEach(() => {
    // Create minimal Express app with just the endpoint
    app = express();
    app.use(express.json());
    
    // Mock structuralizer
    mockStructuralizer = {
      chemicals: jest.fn().mockResolvedValue({
        object: "water",
        chemicals: [
          { name: "Water", smiles: "O" }
        ]
      })
    };
    
    // Replicate the actual route logic
    app.post('/api/structuralize', async (req, res) => {
      const { text, lookupMode = 'GPT-5' } = req.body;

      if (!text || typeof text !== 'string') {
        return res.status(400).json({
          error: 'Missing or invalid "text" parameter'
        });
      }

      if (typeof lookupMode !== 'string') {
        return res.status(400).json({
          error: 'Invalid "lookupMode" parameter'
        });
      }

      try {
        const result = await mockStructuralizer.chemicals({
          object: text.trim(),
          lookupMode
        });
        res.json(result);
      } catch (error) {
        res.status(500).json({ error: error.message });
      }
    });
  });

  describe('Missing/Invalid text parameter', () => {
    test('should return 400 when text is missing', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({});

      expect(response.status).toBe(400);
      expect(response.body.error).toContain('text');
    });

    test('should return 400 when text is empty string', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: '' });

      expect(response.status).toBe(400);
      expect(response.body.error).toContain('text');
    });

    test('should return 400 when text is null', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: null });

      expect(response.status).toBe(400);
      expect(response.body.error).toContain('text');
    });

    test('should return 400 when text is a number', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: 123 });

      expect(response.status).toBe(400);
      expect(response.body.error).toContain('text');
    });

    test('should return 400 when text is an array', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: ['water'] });

      expect(response.status).toBe(400);
      expect(response.body.error).toContain('text');
    });

    test('should return 400 when text is an object', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: { value: 'water' } });

      expect(response.status).toBe(400);
      expect(response.body.error).toContain('text');
    });
  });

  describe('Invalid lookupMode parameter', () => {
    test('should return 400 when lookupMode is a number', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: 'water', lookupMode: 123 });

      expect(response.status).toBe(400);
      expect(response.body.error).toContain('lookupMode');
    });

    test('should return 400 when lookupMode is an array', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: 'water', lookupMode: ['PubChem'] });

      expect(response.status).toBe(400);
      expect(response.body.error).toContain('lookupMode');
    });
  });

  describe('Valid requests', () => {
    test('should return 200 with valid text', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: 'water' });

      expect(response.status).toBe(200);
      expect(response.body).toHaveProperty('object');
      expect(response.body).toHaveProperty('chemicals');
    });

    test('should accept valid lookupMode string', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: 'coffee', lookupMode: 'PubChem' });

      expect(response.status).toBe(200);
      expect(mockStructuralizer.chemicals).toHaveBeenCalledWith({
        object: 'coffee',
        lookupMode: 'PubChem'
      });
    });

    test('should use default lookupMode when not provided', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: 'salt' });

      expect(response.status).toBe(200);
      expect(mockStructuralizer.chemicals).toHaveBeenCalledWith({
        object: 'salt',
        lookupMode: 'GPT-5'
      });
    });

    test('should trim whitespace from text', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: '  ethanol  ' });

      expect(response.status).toBe(200);
      expect(mockStructuralizer.chemicals).toHaveBeenCalledWith({
        object: 'ethanol',
        lookupMode: 'GPT-5'
      });
    });
  });

  describe('Response format', () => {
    test('should return object and chemicals array', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: 'water' });

      expect(response.status).toBe(200);
      expect(response.body.object).toBe('water');
      expect(Array.isArray(response.body.chemicals)).toBe(true);
      expect(response.body.chemicals.length).toBeGreaterThan(0);
    });

    test('chemicals should have name and smiles properties', async () => {
      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: 'water' });

      expect(response.status).toBe(200);
      response.body.chemicals.forEach(chemical => {
        expect(chemical).toHaveProperty('name');
        expect(chemical).toHaveProperty('smiles');
        expect(typeof chemical.name).toBe('string');
        expect(typeof chemical.smiles).toBe('string');
      });
    });
  });

  describe('Error handling', () => {
    test('should return 500 when structuralizer throws', async () => {
      mockStructuralizer.chemicals.mockRejectedValueOnce(new Error('AI service unavailable'));

      const response = await request(app)
        .post('/api/structuralize')
        .send({ text: 'water' });

      expect(response.status).toBe(500);
      expect(response.body.error).toBe('AI service unavailable');
    });
  });
});
