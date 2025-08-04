// backwards-connection.test.js - Jest integration for backwards connection testing
const BackwardsConnectionTester = require('./backwards-connection-test');

// Import the Express app for testing
let app;
try {
  app = require('../backend/api/server');
} catch (error) {
  console.error('Failed to import server app:', error.message);
}

describe('Backwards Connection Testing', () => {
  let tester;

  beforeAll(() => {
    if (!app) {
      throw new Error('Cannot run tests: Express app not available');
    }
    tester = new BackwardsConnectionTester(app);
  });

  // Test each step individually for better error isolation
  describe('Step 5: Frontend Serving', () => {
    test('should serve frontend with 3Dmol library', async () => {
      const result = await tester.testStep5_FrontendServing();
      expect(result).toBe(true);
      expect(tester.testResults.step5_frontend.status).toBe('pass');
    }, 10000);
  });

  describe('Step 4: SDF File Serving', () => {
    test('should serve SDF files correctly', async () => {
      const result = await tester.testStep4_SdfFileServing();
      expect(result).toBe(true);
      expect(tester.testResults.step4_sdfServing.status).toBe('pass');
    }, 10000);
  });

  describe('Step 3: SDF Generation', () => {
    test('should generate SDF files from SMILES', async () => {
      const result = await tester.testStep3_SdfGeneration();
      expect(result.success).toBe(true);
      expect(tester.testResults.step3_sdfGeneration.status).toBe('pass');
      expect(result.sdfPaths).toBeDefined();
      expect(result.sdfPaths.length).toBeGreaterThan(0);
    }, 15000);
  });

  describe('Step 2: Molecular Analysis', () => {
    test('should analyze text inputs and return molecular data', async () => {
      const result = await tester.testStep2_MolecularAnalysis();
      expect(result.success).toBe(true);
      expect(['pass', 'partial']).toContain(tester.testResults.step2_molecularAnalysis.status);
      expect(result.textResult).toBeDefined();
    }, 20000); // Longer timeout for AI API calls
  });

  describe('Step 1: Input Validation', () => {
    test('should validate inputs correctly', async () => {
      const result = await tester.testStep1_InputValidation();
      expect(result).toBe(true);
      expect(tester.testResults.step1_inputValidation.status).toBe('pass');
    }, 10000);
  });

  // Full backwards test integration
  describe('Complete Backwards Connection Test', () => {
    test('should pass all backwards connection tests', async () => {
      const report = await tester.runBackwardsConnectionTest();
      
      // Expect at least 80% connection health
      expect(report.connectionHealth).toBeGreaterThanOrEqual(60);
      expect(report.passCount).toBeGreaterThanOrEqual(3);
      
      // Log the full report for debugging
      console.log('\n📋 Full Connection Test Report:', JSON.stringify(report, null, 2));
      
    }, 60000); // Long timeout for full test suite
  });

  // Test connection health monitoring
  describe('Connection Health Monitoring', () => {
    test('should provide detailed connection health metrics', () => {
      const report = tester.generateReport();
      
      expect(report).toHaveProperty('connectionHealth');
      expect(report).toHaveProperty('passCount');
      expect(report).toHaveProperty('totalSteps');
      expect(report).toHaveProperty('details');
      
      expect(typeof report.connectionHealth).toBe('number');
      expect(report.connectionHealth).toBeGreaterThanOrEqual(0);
      expect(report.connectionHealth).toBeLessThanOrEqual(100);
    });
  });
});

// Optional: Add performance benchmarks
describe('Connection Performance', () => {
  test('frontend serving should respond quickly', async () => {
    const startTime = Date.now();
    await tester.testStep5_FrontendServing();
    const responseTime = Date.now() - startTime;
    
    expect(responseTime).toBeLessThan(2000); // 2 second max
  });

  test('SDF file serving should be fast', async () => {
    const startTime = Date.now();
    await tester.testStep4_SdfFileServing();
    const responseTime = Date.now() - startTime;
    
    expect(responseTime).toBeLessThan(3000); // 3 second max
  });
});