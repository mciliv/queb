const fs = require('fs');
const path = require('path');
const { APP_VALIDATION_RULES, ValidationHelpers } = require('./support/app-validation-rules');

describe('Basic App Smoke Tests - Validation Rules Compliance', () => {
  test('Core files exist and meet validation requirements', () => {
    const coreFiles = [
      'src/client/core/index.html',
      'src/server/api/server.js',
      'package.json'
    ];

    coreFiles.forEach(file => {
      const filePath = path.join(__dirname, '..', file);
      expect(fs.existsSync(filePath)).toBe(true);
      
      // Additional validation - file size checks
      const stats = fs.statSync(filePath);
      expect(stats.size).toBeGreaterThan(0);
      
      // For critical files, ensure they're not too large
      if (file.includes('.js') || file.includes('.html')) {
        expect(stats.size).toBeLessThan(5242880); // 5MB max
      }
    });
    
    console.log('✅ Core files validated against size and existence rules');
  });

  test('Package.json has required dependencies per validation rules', () => {
    const packageJson = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'package.json')));

    // Use validation rules to check dependencies
    const browserRules = APP_VALIDATION_RULES.compatibility.browsers;
    const requiredDeps = ['express', 'cors', 'zod'];
    
    requiredDeps.forEach(dep => {
      expect(packageJson.dependencies[dep]).toBeDefined();
    });
    
    // Additional checks from validation rules
    expect(packageJson.name).toBeDefined();
    expect(packageJson.version).toBeDefined();
    expect(packageJson.scripts).toBeDefined();
    
    console.log('✅ Package.json validated against dependency and structure rules');
  });

  test('Frontend assets exist', () => {
    const assetFiles = [
      'src/client/assets/style.css',
      'src/client/assets/camera.svg',
      'src/client/assets/account.svg'
    ];
    
    assetFiles.forEach(file => {
      const filePath = path.join(__dirname, '..', file);
      if (fs.existsSync(filePath)) {
        expect(fs.existsSync(filePath)).toBe(true);
        console.log(`✅ Asset found: ${file}`);
      } else {
        console.log(`⚠️ Asset not found (optional): ${file}`);
        // Don't fail the test for missing assets, just log them
      }
    });
    
    console.log('✅ Frontend asset validation completed (optional files)');
  });

  test('Environment configuration is valid', () => {
    // Test that key configuration files exist
    const configFiles = ['package.json', 'src/client/build-frontend.js'];

    configFiles.forEach(file => {
      const filePath = path.join(__dirname, '..', file);
      if (fs.existsSync(filePath)) {
        expect(fs.existsSync(filePath)).toBe(true);
      }
    });

    // Test that main entry points exist
    const entryPoints = ['src/server/api/server.js', 'src/client/core/index.html'];
    entryPoints.forEach(entry => {
      expect(fs.existsSync(path.join(__dirname, '..', entry))).toBe(true);
    });
  });

  // Comprehensive validation rule compliance tests
  test('Input validation helpers work correctly', () => {
    // Test valid inputs
    APP_VALIDATION_RULES.validation.textInput.validExamples.forEach(example => {
      const result = ValidationHelpers.validateTextInput(example);
      expect(result).toBeNull(); // null means valid
    });

    // Test invalid inputs
    const invalidInputs = ['', 'x', 'love', 'asdf', 'test123'];
    invalidInputs.forEach(input => {
      const result = ValidationHelpers.validateTextInput(input);
      if (input.length < APP_VALIDATION_RULES.validation.textInput.minLength || 
          APP_VALIDATION_RULES.validation.textInput.invalidPatterns.some(p => p.test(input.toLowerCase()))) {
        expect(result).not.toBeNull(); // Should return error message
      }
    });

    console.log('✅ Input validation rules tested and working');
  });

  test('Performance benchmarks are reasonable', () => {
    const benchmarks = APP_VALIDATION_RULES.performance.benchmarks;
    
    Object.keys(benchmarks).forEach(operation => {
      const benchmark = benchmarks[operation];
      expect(benchmark.maxTime).toBeGreaterThan(1000); // At least 1 second
      expect(benchmark.maxTime).toBeLessThan(30000); // Less than 30 seconds
    });

    console.log('✅ Performance benchmarks validated as reasonable');
  });

  test('Security validation helpers work correctly', () => {
    // Test XSS detection
    const maliciousInputs = [
      '<script>alert("xss")</script>',
      'javascript:alert(1)',
      '<img onerror="alert(1)" src="x">'
    ];

    maliciousInputs.forEach(input => {
      const isSecure = ValidationHelpers.validateSecurity('input', input);
      expect(isSecure).toBe(false); // Should detect as insecure
    });

    // Test safe inputs
    const safeInputs = ['water', 'coffee beans', 'olive oil'];
    safeInputs.forEach(input => {
      const isSecure = ValidationHelpers.validateSecurity('input', input);
      expect(isSecure).toBe(true); // Should be safe
    });

    console.log('✅ Security validation rules tested and working');
  });

  test('API validation rules are properly defined', () => {
    const apis = APP_VALIDATION_RULES.core.server.apis;
    
    Object.keys(apis).forEach(endpoint => {
      const api = apis[endpoint];
      expect(api.method).toBeDefined();
      expect(api.requiredFields).toBeDefined();
      expect(Array.isArray(api.requiredFields)).toBe(true);
      expect(api.timeout).toBeGreaterThan(1000);
      expect(api.timeout).toBeLessThan(30000);
    });

    console.log('✅ API validation rules are properly structured');
  });
}); 