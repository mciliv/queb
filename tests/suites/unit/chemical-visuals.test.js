/**
 * This file previously contained an unscoped snippet referencing `viewer`,
 * which caused unit test runs to fail with ReferenceError.
 *
 * The actual 3D viewer (3Dmol.js) is browser-only and belongs in visual/e2e tests.
 */

describe('Chemical visuals (placeholder)', () => {
  test('placeholder to prevent empty suite error', () => {
    expect(true).toBe(true);
  });
});