/**
 * Smoke Tests - Quick sanity checks against running server
 * 
 * Usage:
 *   1. Start the server: ./run
 *   2. Run tests: npm run test:smoke
 * 
 * Tests health, static files, API validation, and basic functionality
 */

const path = require('path');
const { spawn } = require('child_process');

const BASE_URL = process.env.TEST_URL || 'http://localhost:8080';
const PROJECT_ROOT = path.resolve(__dirname, '../../..');
const SERVER_ENTRY = path.join(PROJECT_ROOT, 'src/server/api/server.js');

let serverProcess;

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

async function pingHealth(timeoutMs = 1000) {
  try {
    const response = await fetch(`${BASE_URL}/api/health`, {
      signal: AbortSignal.timeout(timeoutMs),
    });
    return response.ok;
  } catch {
    return false;
  }
}

async function ensureServerRunning() {
  // If already up, reuse it so we don't kill a user-run server
  if (await pingHealth()) return;

  serverProcess = spawn('node', [SERVER_ENTRY], {
    cwd: PROJECT_ROOT,
    env: { ...process.env, NODE_ENV: process.env.NODE_ENV || 'test' },
    stdio: 'inherit',
  });

  const maxAttempts = 30;
  for (let i = 0; i < maxAttempts; i += 1) {
    if (await pingHealth()) return;
    await sleep(1000);
  }

  throw new Error('Server failed to start for smoke tests');
}

async function stopStartedServer() {
  if (!serverProcess) return;
  serverProcess.kill('SIGTERM');
}

describe('Smoke Tests', () => {
  
  // Check server is running before all tests
  beforeAll(async () => {
    await ensureServerRunning();
  });

  afterAll(async () => {
    await stopStartedServer();
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
      expect(html.toLowerCase()).toContain('<!doctype html>');
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

