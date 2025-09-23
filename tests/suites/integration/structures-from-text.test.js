// E2E test for the unified text → structures flow
// Drives the new /structures-from-text endpoint end-to-end

const request = require('supertest');
const fs = require('fs');
const path = require('path');

// Use the Express app directly (local server impl)
const app = require('../../../src/server/api/server');

describe('End-to-End: /structures-from-text', () => {
  beforeAll(() => {
    process.env.NODE_ENV = 'test';
  });

  const cases = [
    { label: 'black coffee', requiredOneOf: ['caffeine', 'water'], min: 2 },
    { label: 'water', requiredOneOf: ['water'], min: 1 },
  ];

  cases.forEach(({ label, requiredOneOf, min }) => {
    test(`should return structures for "${label}"`, async () => {
      const res = await request(app)
        .post('/structures-from-text')
        .send({ object: label })
        .expect(200);

      // Response shape (no wrapper)
      expect(res.body).toBeDefined();
      const chemicals = res.body.chemicals || res.body.molecules || [];
      expect(Array.isArray(chemicals)).toBe(true);
      expect(chemicals.length).toBeGreaterThanOrEqual(min);

      // Check at least one required name appears (case-insensitive substring)
      const names = chemicals.map(c => (c.name || '').toLowerCase());
      const hit = requiredOneOf.some(req => names.some(n => n.includes(req)));
      expect(hit).toBe(true);

      // Validate SDF accessibility for first up to 3 with sdfPath
      const withSdf = chemicals.filter(c => !!c.sdfPath).slice(0, 3);
      for (const item of withSdf) {
        const httpResp = await request(app).get(item.sdfPath);
        expect(httpResp.status).toBe(200);
        const text = httpResp.text || '';
        expect(text.length).toBeGreaterThan(0);
        expect(text.includes('$$$$')).toBe(true);

        // Also confirm file exists on disk (served from test/sdf_files)
        const fileName = item.sdfPath.replace('/sdf_files/', '');
        const full = path.join(__dirname, '../../../test/sdf_files', fileName);
        expect(fs.existsSync(full)).toBe(true);
      }
    }, 30000);
  });
});


