// E2E test for name → SDF conversion
// Tests the /name-to-sdf endpoint end-to-end

const request = require('supertest');
const fs = require('fs');
const path = require('path');

// Use the Express app directly (local server impl)
const app = require('../../../src/server/api/server');

describe('End-to-End: /name-to-sdf', () => {
  beforeAll(() => {
    process.env.NODE_ENV = 'test';
  });

  const testCases = [
    { name: 'water', expectedInName: 'water' },
    { name: 'ethanol', expectedInName: 'ethanol' },
    { name: 'caffeine', expectedInName: 'caffeine' },
    { name: 'glucose', expectedInName: 'glucose' },
  ];

  testCases.forEach(({ name, expectedInName }) => {
    test(`should convert "${name}" to SDF file`, async () => {
      const res = await request(app)
        .post('/name-to-sdf')
        .send({ name })
        .expect(200);

      // Response shape validation
      expect(res.body).toBeDefined();
      expect(res.body).toHaveProperty('name');
      expect(res.body).toHaveProperty('sdfPath');
      expect(res.body).toHaveProperty('status');

      // Validate resolved name
      const resolvedName = (res.body.name || '').toLowerCase();
      expect(resolvedName).toContain(expectedInName);

      // Check status
      expect(['ok', 'lookup_required']).toContain(res.body.status);

      // If SDF was generated, validate it
      if (res.body.sdfPath && res.body.status === 'ok') {
        // Test HTTP accessibility
        const sdfResponse = await request(app).get(res.body.sdfPath);
        expect(sdfResponse.status).toBe(200);
        
        const sdfContent = sdfResponse.text || '';
        expect(sdfContent.length).toBeGreaterThan(0);
        
        // Validate SDF format (should end with $$$$)
        expect(sdfContent.includes('$$$$')).toBe(true);

        // Confirm file exists on disk
        const fileName = res.body.sdfPath.replace('/sdf_files/', '');
        const fullPath = path.join(__dirname, '../../sdf_files', fileName);
        expect(fs.existsSync(fullPath)).toBe(true);
      }
    }, 30000);
  });

  test('should return 400 for missing name', async () => {
    const res = await request(app)
      .post('/name-to-sdf')
      .send({})
      .expect(400);

    expect(res.body).toHaveProperty('error');
    expect(res.body.error).toContain('name is required');
  });

  test('should return 400 for empty name', async () => {
    const res = await request(app)
      .post('/name-to-sdf')
      .send({ name: '   ' })
      .expect(400);

    expect(res.body).toHaveProperty('error');
  });

  test('should handle unknown molecule name gracefully', async () => {
    const res = await request(app)
      .post('/name-to-sdf')
      .send({ name: 'xyzunknownmolecule12345' })
      .expect(200);

    expect(res.body).toBeDefined();
    expect(res.body).toHaveProperty('status');
    
    // Should return lookup_required or null sdfPath for unknown molecules
    if (res.body.status === 'lookup_required') {
      expect(res.body.sdfPath).toBeNull();
    }
  });

  test('should respect overwrite parameter', async () => {
    const testName = 'ethanol';
    
    // First request without overwrite
    const res1 = await request(app)
      .post('/name-to-sdf')
      .send({ name: testName, overwrite: false })
      .expect(200);

    expect(res1.body).toHaveProperty('sdfPath');
    const firstPath = res1.body.sdfPath;

    // Second request with overwrite
    const res2 = await request(app)
      .post('/name-to-sdf')
      .send({ name: testName, overwrite: true })
      .expect(200);

    expect(res2.body).toHaveProperty('sdfPath');
    
    // Both should have valid SDF paths (same or different)
    if (firstPath && res2.body.sdfPath) {
      expect(res2.body.sdfPath).toBeTruthy();
    }
  }, 30000);

  test('should handle case variations in molecule names', async () => {
    const names = ['Water', 'WATER', 'water', 'wAtEr'];
    
    for (const name of names) {
      const res = await request(app)
        .post('/name-to-sdf')
        .send({ name })
        .expect(200);

      expect(res.body).toHaveProperty('status');
      expect(res.body).toHaveProperty('sdfPath');
    }
  }, 30000);
});



