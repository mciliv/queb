// Simple smoke tests - validate basic app functionality
const fs = require('fs');
const path = require('path');

describe('Basic App Smoke Tests', () => {
  test('Core files exist', () => {
    const coreFiles = [
      'frontend/core/index.html',
      'frontend/core/app.js', 
      'backend/api/server.js',
      'package.json'
    ];
    
    coreFiles.forEach(file => {
      expect(fs.existsSync(path.join(__dirname, '..', file))).toBe(true);
    });
  });

  test('Package.json has required dependencies', () => {
    const packageJson = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'package.json')));
    
    const requiredDeps = ['express', 'cors', 'openai'];
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
    // Test that root-level scripts are set up
    const rootScripts = ['dev', 'tests', 'ship', 'server', 'debug', 'cleanup'];
    
    rootScripts.forEach(script => {
      expect(fs.existsSync(path.join(__dirname, '..', script))).toBe(true);
    });
    
    // Test that infrastructure scripts still exist
    expect(fs.existsSync(path.join(__dirname, '..', 'infrastructure', 'scripts'))).toBe(true);
  });
}); 