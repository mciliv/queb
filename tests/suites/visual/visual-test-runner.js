/**
 * Visual Test Runner
 * Runs automated visual tests for molecular visualization
 * Usage: node tests/suites/visual/visual-test-runner.js
 */

const { processSmilesToViewers } = require('../../../src/client/utils/visual-test-helpers.js');
const { PRESET_VISUAL_TESTS, SMILES_NAME_MAP } = require('../../../src/client/constants.js');

/**
 * Simulates the SDF generation API call
 * In a real test, this would make an actual HTTP request
 */
async function generateSDFs(smilesArray) {
  // Mock implementation - replace with actual API call in integration tests
  const sanitizeSmiles = (s) => s.replace(/[^a-zA-Z0-9]/g, ch => ch === '=' ? '__' : '_');
  
  return {
    sdfPaths: smilesArray.map(sm => `/sdf_files/${sanitizeSmiles(sm)}.sdf`)
  };
}

/**
 * Runs a single visual test
 */
async function runVisualTest(test) {
  const { label, smilesList } = test;
  
  if (!smilesList || smilesList.length === 0) {
    console.warn(`[SKIP] Test "${label}": No SMILES provided`);
    return { success: false, reason: 'no_smiles' };
  }
  
  try {
    // Generate SDFs
    const sdfResult = await generateSDFs(smilesList);
    const returnedPaths = Array.isArray(sdfResult?.sdfPaths) ? sdfResult.sdfPaths : [];
    
    // Process into viewer objects
    const viewers = processSmilesToViewers(smilesList, returnedPaths, SMILES_NAME_MAP);
    
    // Validate results
    const valid = viewers.every(v => 
      v.smiles && 
      v.name && 
      v.sdfData && 
      v.sdfData.startsWith('file://')
    );
    
    if (!valid) {
      console.error(`[FAIL] Test "${label}": Invalid viewer objects`);
      return { success: false, reason: 'invalid_viewers', viewers };
    }
    
    console.log(`[PASS] Test "${label}": Generated ${viewers.length} viewers`);
    return { success: true, viewers };
    
  } catch (error) {
    console.error(`[ERROR] Test "${label}":`, error.message);
    return { success: false, reason: 'exception', error: error.message };
  }
}

/**
 * Runs all preset visual tests
 */
async function runAllVisualTests() {
  console.log('='.repeat(60));
  console.log('Visual Test Runner');
  console.log('='.repeat(60));
  console.log(`Running ${PRESET_VISUAL_TESTS.length} preset tests...\n`);
  
  const results = [];
  
  for (const test of PRESET_VISUAL_TESTS) {
    const result = await runVisualTest(test);
    results.push({ test: test.label, ...result });
  }
  
  // Summary
  console.log('\n' + '='.repeat(60));
  console.log('Test Summary');
  console.log('='.repeat(60));
  
  const passed = results.filter(r => r.success).length;
  const failed = results.length - passed;
  
  console.log(`Total: ${results.length}`);
  console.log(`Passed: ${passed}`);
  console.log(`Failed: ${failed}`);
  
  if (failed > 0) {
    console.log('\nFailed tests:');
    results.filter(r => !r.success).forEach(r => {
      console.log(`  - ${r.test}: ${r.reason}`);
    });
  }
  
  return { passed, failed, results };
}

// Run if called directly
if (require.main === module) {
  runAllVisualTests()
    .then(({ passed, failed }) => {
      process.exit(failed > 0 ? 1 : 0);
    })
    .catch(error => {
      console.error('Fatal error:', error);
      process.exit(1);
    });
}

module.exports = { runVisualTest, runAllVisualTests };

