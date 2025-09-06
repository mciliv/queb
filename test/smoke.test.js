// Simple smoke tests - validate basic app functionality
const fs = require('fs');
const path = require('path');

describe('Basic App Smoke Tests', () => {
  test('Core files exist', () => {
    const coreFiles = [
      'frontend/core/index.html',
      'backend/api/server.js',
      'package.json'
    ];

    coreFiles.forEach(file => {
      expect(fs.existsSync(path.join(__dirname, '..', file))).toBe(true);
    });
  });

  test('Package.json has required dependencies', () => {
    const packageJson = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'package.json')));

    const requiredDeps = ['express', 'cors'];
    requiredDeps.forEach(dep => {
      expect(packageJson.dependencies[dep]).toBeDefined();
    });
  });

  test('Frontend assets exist', () => {
    const assetFiles = [
      'frontend/assets/style.css',
      'frontend/assets/camera.svg',
      'frontend/assets/account.svg'
    ];
    
    assetFiles.forEach(file => {
      expect(fs.existsSync(path.join(__dirname, '..', file))).toBe(true);
    });
  });

  test('Environment configuration is valid', () => {
    // Test that key configuration files exist
    const configFiles = ['package.json', 'frontend/build-frontend.js'];

    configFiles.forEach(file => {
      expect(fs.existsSync(path.join(__dirname, '..', file))).toBe(true);
    });

    // Test that main entry point exists
    expect(fs.existsSync(path.join(__dirname, '..', 'index.js'))).toBe(true);
  });
}); 