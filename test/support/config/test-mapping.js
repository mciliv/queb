/**
 * File-Based Test Mapping Configuration
 * 
 * This file defines which tests should run when specific files change.
 * Tests are organized by dependency layers to solve dependency issues.
 */

const path = require('path');

/**
 * Test Categories by Dependency Level
 * Lower levels can depend on higher levels, but not vice versa
 */
const TEST_LAYERS = {
  // Level 1: Pure utility functions, no dependencies
  UTILS: 'utils',
  
  // Level 2: Individual modules/services with minimal dependencies
  UNIT: 'unit',
  
  // Level 3: API endpoints and service integration
  API: 'api',
  
  // Level 4: Component integration (frontend + backend)
  INTEGRATION: 'integration',
  
  // Level 5: End-to-end application tests
  E2E: 'e2e'
};

/**
 * File patterns to test mappings
 * Each entry maps file patterns to test categories and specific test files
 */
const FILE_TEST_MAPPING = {
  // Backend Services
  'backend/services/Structuralizer.js': {
    direct: ['test/unit/unit.test.js'],
    layers: [TEST_LAYERS.UNIT, TEST_LAYERS.API, TEST_LAYERS.INTEGRATION],
    related: [
      'test/unit/prompt-accuracy.test.js',
      'test/integration/molecular-accuracy.test.js',
      'backend/test-structuralizer.js'
    ]
  },

  // Test Configuration
  'backend/test-config-example.js': {
    direct: ['backend/test-structuralizer.js', 'test/unit/unit.test.js'],
    layers: [TEST_LAYERS.UNIT],
    related: []
  },

  'backend/test-structuralizer.js': {
    direct: ['test/unit/unit.test.js'],
    layers: [TEST_LAYERS.UNIT],
    related: ['backend/test-config-example.js']
  },

  'backend/services/molecular-processor.js': {
    direct: ['test/unit/unit.test.js'],
    layers: [TEST_LAYERS.UNIT, TEST_LAYERS.API],
    related: []
  },

  'backend/services/user-service.js': {
    direct: ['test/unit/unit.test.js'],
    layers: [TEST_LAYERS.UNIT, TEST_LAYERS.API],
    related: []
  },

  // Backend API
  'backend/api/server.js': {
    direct: ['test/unit/unit.test.js', 'test/integration/server.test.js'],
    layers: [TEST_LAYERS.API, TEST_LAYERS.INTEGRATION],
    related: [
      'test/integration/smoke.test.js'
    ]
  },

  // Backend Schemas
  'backend/schemas/schemas.js': {
    direct: ['test/unit/unit.test.js'],
    layers: [TEST_LAYERS.UNIT, TEST_LAYERS.API],
    related: []
  },

  // Backend Prompts
  'backend/prompts/**.js': {
    direct: ['test/unit/unit.test.js'],
    layers: [TEST_LAYERS.UNIT],
    related: [
      'test/unit/prompt-accuracy.test.js',
      'test/integration/molecular-accuracy.test.js'
    ]
  },

  // Frontend Core
  'frontend/core/App.jsx': {
    direct: ['test/unit/manual.test.js'],
    layers: [TEST_LAYERS.UNIT, TEST_LAYERS.INTEGRATION, TEST_LAYERS.E2E],
    related: [
      'test/integration/integration.test.js'
    ]
  },

  'frontend/core/index.html': {
    direct: [],
    layers: [TEST_LAYERS.INTEGRATION, TEST_LAYERS.E2E],
    related: [
      'test/integration/integration.test.js',
      'test/integration/smoke.test.js'
    ]
  },

  // Frontend Components
  'frontend/components/camera.js': {
    direct: ['test/unit/camera.test.js'],
    layers: [TEST_LAYERS.UNIT, TEST_LAYERS.INTEGRATION],
    related: [
      'test/integration/integration.test.js'
    ]
  },

  'frontend/components/camera-handler.js': {
    direct: ['test/unit/camera-handler.test.js'],
    layers: [TEST_LAYERS.UNIT, TEST_LAYERS.INTEGRATION],
    related: [
      'test/unit/camera.test.js',
      'test/integration/integration.test.js'
    ]
  },

  'frontend/components/payment.js': {
    direct: [],
    layers: [TEST_LAYERS.UNIT, TEST_LAYERS.INTEGRATION],
    related: [
      'test/integration/integration.test.js'
    ]
  },

  'frontend/components/error-handler.js': {
    direct: [],
    layers: [TEST_LAYERS.UNIT, TEST_LAYERS.INTEGRATION],
    related: [
      'test/integration/integration.test.js'
    ]
  },

  // Frontend Assets
  'frontend/assets/style.css': {
    direct: [],
    layers: [TEST_LAYERS.INTEGRATION],
    related: [
      'test/integration/smoke.test.js'
    ]
  },

  // Chemistry Processors (Python removed)

  // Visual/E2E Tests
  'test/suites/e2e/automated-visual-tests.test.js': {
    direct: [],
    layers: [TEST_LAYERS.E2E],
    related: [
      'backend/test-config-example.js',
      'test/integration/smoke.test.js'
    ]
  },

  // Infrastructure
  'package.json': {
    direct: [],
    layers: [TEST_LAYERS.UNIT, TEST_LAYERS.API, TEST_LAYERS.INTEGRATION, TEST_LAYERS.E2E],
    related: [
      'test/integration/smoke.test.js'
    ]
  }
};

/**
 * Test suites organized by dependency layers
 */
const LAYER_TEST_SUITES = {
  [TEST_LAYERS.UTILS]: [
    // Pure utility tests would go here
  ],

  [TEST_LAYERS.UNIT]: [
    'test/unit/unit.test.js',
    'test/unit/camera.test.js', 
    'test/unit/camera-handler.test.js',
    'test/unit/manual.test.js',
    'test/unit/prompt-accuracy.test.js',
    'test/unit/test_*.py'
  ],

  [TEST_LAYERS.API]: [
    'test/integration/server.test.js'
  ],

  [TEST_LAYERS.INTEGRATION]: [
    'test/integration/integration.test.js',
    'test/integration/molecular-accuracy.test.js',
    'test/integration/smoke.test.js',
    'test/integration/system.test.js'
  ],

  [TEST_LAYERS.E2E]: [
    'test/suites/e2e/automated-visual-tests.test.js',
    'test/integration/deployment.test.js'
  ]
};

/**
 * Get tests to run for a changed file
 */
function getTestsForFile(filePath) {
  const normalizedPath = path.normalize(filePath);
  let testsToRun = new Set();

  // Check direct mappings
  for (const [pattern, config] of Object.entries(FILE_TEST_MAPPING)) {
    if (matchesPattern(normalizedPath, pattern)) {
      // Add direct tests
      config.direct.forEach(test => testsToRun.add(test));
      
      // Add related tests
      config.related.forEach(test => testsToRun.add(test));
      
      // Add layer tests
      config.layers.forEach(layer => {
        LAYER_TEST_SUITES[layer].forEach(test => testsToRun.add(test));
      });
    }
  }

  return Array.from(testsToRun);
}

/**
 * Simple pattern matching (supports ** for directories)
 */
function matchesPattern(filePath, pattern) {
  const regex = pattern
    .replace(/\*\*/g, '.*')
    .replace(/\*/g, '[^/]*')
    .replace(/\./g, '\\.');
  
  return new RegExp(`^${regex}$`).test(filePath);
}

/**
 * Get all tests for a specific layer and below
 */
function getTestsForLayer(layer) {
  const layerOrder = Object.values(TEST_LAYERS);
  const currentIndex = layerOrder.indexOf(layer);
  
  let tests = new Set();
  
  for (let i = 0; i <= currentIndex; i++) {
    const currentLayer = layerOrder[i];
    LAYER_TEST_SUITES[currentLayer].forEach(test => tests.add(test));
  }
  
  return Array.from(tests);
}

module.exports = {
  TEST_LAYERS,
  FILE_TEST_MAPPING,
  LAYER_TEST_SUITES,
  getTestsForFile,
  getTestsForLayer,
  matchesPattern
};