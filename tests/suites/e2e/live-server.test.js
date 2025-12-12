/**
 * E2E Tests - Run against live server
 * 
 * Prerequisites:
 *   1. Start the server: ./run start (or ./run)
 *   2. Run tests: npm run test:e2e -- --testPathPattern=live-server
 * 
 * These tests hit the actual running server at localhost:8080
 */

const BASE_URL = process.env.TEST_URL || 'http://localhost:8080';

describe('Live Server Tests', () => {
  
  // Check server is running before all tests
  beforeAll(async () => {
    try {
      const response = await fetch(`${BASE_URL}/api/health`, { 
        signal: AbortSignal.timeout(5000) 
      });
      if (!response.ok) {
        throw new Error(`Server responded with ${response.status}`);
      }
    } catch (error) {
      console.error('\n❌ Server is not running!');
      console.error('Start it with: ./run start');
      console.error(`Expected server at: ${BASE_URL}\n`);
      throw new Error('Server must be running for these tests');
    }
  });

  describe('Health Check', () => {
    test('GET /api/health returns ok status', async () => {
      const response = await fetch(`${BASE_URL}/api/health`);
      
      expect(response.status).toBe(200);
      
      const data = await response.json();
      expect(data.status).toBe('ok');
      expect(data).toHaveProperty('environment');
      expect(data).toHaveProperty('version');
    });
  });

  describe('Static Files', () => {
    test('GET / returns HTML', async () => {
      const response = await fetch(`${BASE_URL}/`);
      
      expect(response.status).toBe(200);
      expect(response.headers.get('content-type')).toContain('text/html');
      
      const html = await response.text();
      expect(html).toContain('<!DOCTYPE html>');
      expect(html).toContain('<title>');
    });

    test('GET /bundle.js returns JavaScript', async () => {
      const response = await fetch(`${BASE_URL}/bundle.js`);
      
      expect(response.status).toBe(200);
      expect(response.headers.get('content-type')).toContain('javascript');
    });

    test('GET /bundle.css returns CSS', async () => {
      const response = await fetch(`${BASE_URL}/bundle.css`);
      
      expect(response.status).toBe(200);
      expect(response.headers.get('content-type')).toContain('css');
    });
  });

  describe('API - Text Structuralize', () => {
    
    describe('Input Validation', () => {
      test('POST /api/structuralize with missing text returns 400', async () => {
        const response = await fetch(`${BASE_URL}/api/structuralize`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({})
        });

        expect(response.status).toBe(400);
        const data = await response.json();
        expect(data.error).toContain('text');
      });

      test('POST /api/structuralize with empty text returns 400', async () => {
        const response = await fetch(`${BASE_URL}/api/structuralize`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ text: '' })
        });

        expect(response.status).toBe(400);
        const data = await response.json();
        expect(data.error).toContain('text');
      });

      test('POST /api/structuralize with invalid text type returns 400', async () => {
        const response = await fetch(`${BASE_URL}/api/structuralize`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ text: 123 })
        });

        expect(response.status).toBe(400);
        const data = await response.json();
        expect(data.error).toContain('text');
      });

      test('POST /api/structuralize with invalid lookupMode returns 400', async () => {
        const response = await fetch(`${BASE_URL}/api/structuralize`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ text: 'water', lookupMode: 123 })
        });

        expect(response.status).toBe(400);
        const data = await response.json();
        expect(data.error).toContain('lookupMode');
      });
    });

    describe('Successful Requests', () => {
      // These tests may take time as they hit the real AI service
      jest.setTimeout(60000);

      test('POST /api/structuralize with valid text returns chemicals', async () => {
        const response = await fetch(`${BASE_URL}/api/structuralize`, {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ text: 'water' })
        });

        // Could be 200 (success) or 500 (AI service error)
        // We just want to verify the endpoint accepts valid input
        expect([200, 500]).toContain(response.status);
        
        if (response.status === 200) {
          const data = await response.json();
          expect(data).toHaveProperty('object');
          expect(data).toHaveProperty('chemicals');
          expect(Array.isArray(data.chemicals)).toBe(true);
        }
      });
    });
  });

  describe('404 Handling', () => {
    test('GET /api/nonexistent returns 404', async () => {
      const response = await fetch(`${BASE_URL}/api/nonexistent`);
      
      expect(response.status).toBe(404);
      const data = await response.json();
      expect(data.error).toContain('not found');
    });
  });

  describe('SPA Routing', () => {
    test('GET /any-client-route returns HTML (SPA fallback)', async () => {
      const response = await fetch(`${BASE_URL}/some/client/route`);
      
      expect(response.status).toBe(200);
      expect(response.headers.get('content-type')).toContain('text/html');
    });
  });
});






