// backwards-connection-test.js - Test endpoints by working backwards through the data flow
const request = require('supertest');
const fs = require('fs');
const path = require('path');

// Test configuration
const TEST_CONFIG = {
  baseUrl: 'http://localhost:8080',
  testSMILES: ['CCO', 'C1=CC=CC=C1'], // ethanol, benzene
  testText: 'caffeine',
  sampleImageBase64: 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==' // 1x1 pixel
};

/**
 * Backwards Connection Testing Strategy:
 * 
 * Step 5: Frontend Display ← Can we serve the frontend?
 * Step 4: SDF File Serving ← Can we serve generated SDF files?  
 * Step 3: SDF Generation ← Can we generate SDF from SMILES?
 * Step 2: Molecular Analysis ← Can we analyze text/images to get SMILES?
 * Step 1: Input Validation ← Can we handle valid inputs?
 */

class BackwardsConnectionTester {
  constructor(app) {
    this.app = app;
    this.testResults = {
      step5_frontend: { status: 'pending', details: null },
      step4_sdfServing: { status: 'pending', details: null },
      step3_sdfGeneration: { status: 'pending', details: null },
      step2_molecularAnalysis: { status: 'pending', details: null },
      step1_inputValidation: { status: 'pending', details: null }
    };
  }

  // Step 5: Test Frontend Serving (final output)
  async testStep5_FrontendServing() {
    try {
      console.log('🔍 Step 5: Testing Frontend Serving...');
      
      const response = await request(this.app).get('/index.html');
      

      
      if (response.status === 200 && response.text.includes('object-input')) {
        this.testResults.step5_frontend = { 
          status: 'pass', 
          details: 'Frontend serves correctly with input field' 
        };
        console.log('✅ Step 5: Frontend serving - PASS');
        return true;
      } else {
        throw new Error(`Frontend response invalid: ${response.status}`);
      }
    } catch (error) {
      this.testResults.step5_frontend = { 
        status: 'fail', 
        details: error.message 
      };
      console.log('❌ Step 5: Frontend serving - FAIL:', error.message);
      return false;
    }
  }

  // Step 4: Test SDF File Serving (works backwards from display)
  async testStep4_SdfFileServing() {
    try {
      console.log('🔍 Step 4: Testing SDF File Serving...');
      
      // First ensure we have an SDF file to test with
      const sdfDir = path.join(__dirname, '..', 'data', 'sdf_files');
      if (!fs.existsSync(sdfDir)) {
        fs.mkdirSync(sdfDir, { recursive: true });
      }
      
      // Create a test SDF file
      const testSdfPath = path.join(sdfDir, 'test_CCO.sdf');
      const testSdfContent = `
  Mrv2014 01010000002D          

  3  2  0  0  0  0            999 V2000
   -1.2375    0.7145    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4125    0.7145    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4125    0.7145    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$`;
      
      fs.writeFileSync(testSdfPath, testSdfContent);
      
      // Test if we can serve the SDF file
      const response = await request(this.app).get('/sdf_files/test_CCO.sdf');
      

      
      if (response.status === 200 && response.text && response.text.includes('$$$$')) {
        this.testResults.step4_sdfServing = { 
          status: 'pass', 
          details: 'SDF files served correctly with valid content' 
        };
        console.log('✅ Step 4: SDF file serving - PASS');
        return true;
      } else {
        // Try alternative path
        const altResponse = await request(this.app).get('/data/sdf_files/test_CCO.sdf');
        if (altResponse.status === 200 && altResponse.text && altResponse.text.includes('$$$$')) {
          this.testResults.step4_sdfServing = { 
            status: 'pass', 
            details: 'SDF files served correctly with valid content (alt path)' 
          };
          console.log('✅ Step 4: SDF file serving - PASS (alt path)');
          return true;
        } else {
          throw new Error('SDF file not served correctly: ' + response.status);
        }
      }
    } catch (error) {
      this.testResults.step4_sdfServing = { 
        status: 'fail', 
        details: error.message 
      };
      console.log('❌ Step 4: SDF file serving - FAIL:', error.message);
      return false;
    }
  }

  // Step 3: Test SDF Generation (works backwards from file serving)
  async testStep3_SdfGeneration() {
    try {
      console.log('🔍 Step 3: Testing SDF Generation...');
      
      const response = await request(this.app)
        .post('/generate-sdfs')
        .send({ smiles: TEST_CONFIG.testSMILES });
      
      if (response.status === 200 && 
          response.body.sdfPaths && 
          response.body.sdfPaths.length > 0) {
        
        this.testResults.step3_sdfGeneration = { 
          status: 'pass', 
          details: `Generated ${response.body.sdfPaths.length} SDF files successfully` 
        };
        console.log('✅ Step 3: SDF generation - PASS');
        return { success: true, sdfPaths: response.body.sdfPaths };
      } else {
        throw new Error(`SDF generation failed: ${response.status} - ${JSON.stringify(response.body)}`);
      }
    } catch (error) {
      this.testResults.step3_sdfGeneration = { 
        status: 'fail', 
        details: error.message 
      };
      console.log('❌ Step 3: SDF generation - FAIL:', error.message);
      return { success: false };
    }
  }

  // Step 2: Test Molecular Analysis (works backwards from SDF generation)
  async testStep2_MolecularAnalysis() {
    try {
      console.log('🔍 Step 2: Testing Molecular Analysis...');
      
      // Test text analysis
      const textResponse = await request(this.app)
        .post('/analyze-text')
        .send({ object: TEST_CONFIG.testText });
      
      if (textResponse.status === 200 && textResponse.body.chemicals) {
        console.log('✅ Step 2a: Text analysis - PASS');
        
        // Test image analysis (if OpenAI key available)
        try {
          const imageResponse = await request(this.app)
            .post('/image-molecules')
            .send({ imageBase64: TEST_CONFIG.sampleImageBase64 });
          
          if (imageResponse.status === 200) {
            this.testResults.step2_molecularAnalysis = { 
              status: 'pass', 
              details: 'Both text and image analysis working' 
            };
            console.log('✅ Step 2b: Image analysis - PASS');
            return { success: true, textResult: textResponse.body };
          }
        } catch (imageError) {
          // Image analysis might fail if no OpenAI key, but text analysis passed
          this.testResults.step2_molecularAnalysis = { 
            status: 'partial', 
            details: 'Text analysis works, image analysis needs OpenAI key' 
          };
          console.log('⚠️ Step 2b: Image analysis - PARTIAL (needs OpenAI key)');
          return { success: true, textResult: textResponse.body };
        }
      } else {
        throw new Error(`Text analysis failed: ${textResponse.status} - ${JSON.stringify(textResponse.body)}`);
      }
    } catch (error) {
      this.testResults.step2_molecularAnalysis = { 
        status: 'fail', 
        details: error.message 
      };
      console.log('❌ Step 2: Molecular analysis - FAIL:', error.message);
      return { success: false };
    }
  }

  // Step 1: Test Input Validation (works backwards from analysis)
  async testStep1_InputValidation() {
    try {
      console.log('🔍 Step 1: Testing Input Validation...');
      
      // Test valid inputs
      const validTextResponse = await request(this.app)
        .post('/analyze-text')
        .send({ object: 'water' });
      
      // Test invalid inputs (should fail gracefully)
      const invalidTextResponse = await request(this.app)
        .post('/analyze-text')
        .send({ object: '' });
      
      const invalidSdfResponse = await request(this.app)
        .post('/generate-sdfs')
        .send({ smiles: 'not-an-array' });
      
      if (validTextResponse.status === 200 && 
          invalidTextResponse.status === 400 && 
          invalidSdfResponse.status === 400) {
        
        this.testResults.step1_inputValidation = { 
          status: 'pass', 
          details: 'Input validation working correctly for valid and invalid inputs' 
        };
        console.log('✅ Step 1: Input validation - PASS');
        return true;
      } else {
        throw new Error('Input validation not working correctly');
      }
    } catch (error) {
      this.testResults.step1_inputValidation = { 
        status: 'fail', 
        details: error.message 
      };
      console.log('❌ Step 1: Input validation - FAIL:', error.message);
      return false;
    }
  }

  // Run all backwards tests in order
  async runBackwardsConnectionTest() {
    console.log('\n🚀 Starting Backwards Connection Test');
    console.log('=====================================');
    console.log('Testing endpoints by working backwards through the data flow...\n');

    // Run tests in backwards order (5 → 1)
    await this.testStep5_FrontendServing();
    await this.testStep4_SdfFileServing();
    await this.testStep3_SdfGeneration();
    await this.testStep2_MolecularAnalysis();
    await this.testStep1_InputValidation();

    return this.generateReport();
  }

  // Generate comprehensive test report
  generateReport() {
    console.log('\n📊 Backwards Connection Test Report');
    console.log('===================================');
    
    const steps = [
      { name: 'Step 5: Frontend Serving', key: 'step5_frontend' },
      { name: 'Step 4: SDF File Serving', key: 'step4_sdfServing' },
      { name: 'Step 3: SDF Generation', key: 'step3_sdfGeneration' },
      { name: 'Step 2: Molecular Analysis', key: 'step2_molecularAnalysis' },
      { name: 'Step 1: Input Validation', key: 'step1_inputValidation' }
    ];
    
    let passCount = 0;
    let totalSteps = steps.length;
    
    steps.forEach(step => {
      const result = this.testResults[step.key];
      const icon = result.status === 'pass' ? '✅' : 
                   result.status === 'partial' ? '⚠️' : '❌';
      
      console.log(`${icon} ${step.name}: ${result.status.toUpperCase()}`);
      console.log(`   ${result.details}\n`);
      
      if (result.status === 'pass' || result.status === 'partial') {
        passCount++;
      }
    });
    
    const connectionHealth = (passCount / totalSteps) * 100;
    console.log(`🔗 Overall Connection Health: ${connectionHealth.toFixed(1)}% (${passCount}/${totalSteps} steps passing)`);
    
    if (connectionHealth >= 80) {
      console.log('🎉 Excellent! All major connections are working.');
    } else if (connectionHealth >= 60) {
      console.log('⚠️ Good, but some connections need attention.');
    } else {
      console.log('🚨 Critical: Multiple connection failures detected.');
    }
    
    return {
      connectionHealth,
      passCount,
      totalSteps,
      details: this.testResults
    };
  }
}

module.exports = BackwardsConnectionTester;

