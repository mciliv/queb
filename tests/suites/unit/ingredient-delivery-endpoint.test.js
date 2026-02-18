/**
 * Unit tests for highest-quality ingredient endpoint validation
 */

const request = require('supertest');
const express = require('express');

jest.mock('../../../src/server/services/ingredient-delivery', () => ({
  fetchHighestQualityIngredients: jest.fn()
}));

const { fetchHighestQualityIngredients } = require('../../../src/server/services/ingredient-delivery');

describe('POST /api/ingredients/highest-quality - Input Validation', () => {
  let app;

  beforeEach(() => {
    app = express();
    app.use(express.json());

    app.post('/api/ingredients/highest-quality', async (req, res) => {
      try {
        const { city, lat, lng, limit } = req.body || {};

        if (city && typeof city !== 'string') {
          return res.status(400).json({ error: 'Invalid "city" parameter' });
        }

        const hasLatLng = lat !== undefined || lng !== undefined;
        if (hasLatLng && (typeof lat !== 'number' || typeof lng !== 'number')) {
          return res.status(400).json({ error: 'Invalid "lat"/"lng" parameters' });
        }

        if (!city && !hasLatLng) {
          return res.status(400).json({ error: 'Provide "city" or "lat"/"lng"' });
        }

        if (limit !== undefined && (!Number.isInteger(limit) || limit <= 0)) {
          return res.status(400).json({ error: 'Invalid "limit" parameter' });
        }

        const results = await fetchHighestQualityIngredients({
          city: city ? city.trim() : undefined,
          lat,
          lng,
          limit
        });

        res.json({ results });
      } catch (error) {
        res.status(500).json({ error: error.message });
      }
    });
  });

  test('returns 400 when missing city and lat/lng', async () => {
    const response = await request(app)
      .post('/api/ingredients/highest-quality')
      .send({});

    expect(response.status).toBe(400);
    expect(response.body.error).toContain('Provide');
  });

  test('returns 400 when city is invalid', async () => {
    const response = await request(app)
      .post('/api/ingredients/highest-quality')
      .send({ city: 123 });

    expect(response.status).toBe(400);
    expect(response.body.error).toContain('city');
  });

  test('returns 400 when lat/lng are invalid', async () => {
    const response = await request(app)
      .post('/api/ingredients/highest-quality')
      .send({ lat: '37.7', lng: -122.4 });

    expect(response.status).toBe(400);
    expect(response.body.error).toContain('lat');
  });

  test('returns 400 when limit is invalid', async () => {
    const response = await request(app)
      .post('/api/ingredients/highest-quality')
      .send({ city: 'Seattle', limit: 0 });

    expect(response.status).toBe(400);
    expect(response.body.error).toContain('limit');
  });

  test('returns 200 with valid city request', async () => {
    fetchHighestQualityIngredients.mockResolvedValueOnce([{ name: 'Line-Caught Salmon' }]);

    const response = await request(app)
      .post('/api/ingredients/highest-quality')
      .send({ city: 'Seattle', limit: 1 });

    expect(response.status).toBe(200);
    expect(response.body.results).toHaveLength(1);
    expect(fetchHighestQualityIngredients).toHaveBeenCalledWith({
      city: 'Seattle',
      lat: undefined,
      lng: undefined,
      limit: 1
    });
  });

  test('returns 200 with valid lat/lng request', async () => {
    fetchHighestQualityIngredients.mockResolvedValueOnce([{ name: 'Heirloom Tomatoes' }]);

    const response = await request(app)
      .post('/api/ingredients/highest-quality')
      .send({ lat: 37.77, lng: -122.41 });

    expect(response.status).toBe(200);
    expect(response.body.results).toHaveLength(1);
    expect(fetchHighestQualityIngredients).toHaveBeenCalledWith({
      city: undefined,
      lat: 37.77,
      lng: -122.41,
      limit: undefined
    });
  });
});
