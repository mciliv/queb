/**
 * Smoke tests - basic functionality verification
 */

describe('Smoke Tests', () => {
  test('server can be imported', () => {
    // Test that the server module can be imported without errors
    expect(() => {
      require('../../server');
    }).not.toThrow();
  });

  test('environment is properly configured', () => {
    // Test that basic environment variables exist
    expect(process.env.NODE_ENV).toBeDefined();
    expect(typeof process.env.PORT).toBe('string');
  });

  test('core services are available', () => {
    // Test that core services can be imported
    expect(() => {
      require('../core/services');
    }).not.toThrow();
  });
});