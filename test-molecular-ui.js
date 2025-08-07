#!/usr/bin/env node

/**
 * Simple molecular visualization UI tester
 * Tests the complete flow: analysis -> SDF generation -> visualization data
 * Run this while the dev server is running
 */

const readline = require('readline');

// Test objects for molecular analysis (matches integration test data)
const testObjects = [
  // Basics - high accuracy expected
  { name: '💧 Water', input: 'water', expected: ['water'], category: 'basics' },
  { name: '🧪 Ethanol', input: 'ethanol', expected: ['ethanol'], category: 'basics' },
  { name: '🧂 Salt', input: 'sodium chloride', expected: ['sodium', 'chloride'], category: 'basics' },
  
  // Beverages - realistic composition
  { name: '🍷 Wine', input: 'red wine', expected: ['ethanol', 'water', 'tartaric acid'], category: 'beverages' },
  { name: '☕ Coffee', input: 'black coffee', expected: ['water', 'caffeine'], category: 'beverages' },
  
  // Biological - complex but realistic
  { name: '🍎 Apple', input: 'fresh apple', expected: ['water', 'fructose', 'glucose'], category: 'biological' },
  
  // Additional useful tests
  { name: '🥬 Kale', input: 'kale', expected: ['water', 'glucose', 'chlorophyll', 'cellulose'], category: 'biological' },
  { name: '🥛 Milk', input: 'milk', expected: ['water', 'lactose', 'proteins', 'fats'], category: 'biological' }
];

async function testMolecularAnalysis(testObject) {
  const baseUrl = 'http://localhost:3000';
  
  try {
    console.log(`\n🧪 Testing: ${testObject.name}`);
    console.log(`   Input: "${testObject.input}"`);
    console.log(`   Category: ${testObject.category || 'general'}`);
    console.log(`   Expected molecules: ${testObject.expected.join(', ')}`);
    
    // Step 1: Analyze the object
    console.log('   📊 Analyzing...');
    const analysisResponse = await fetch(`${baseUrl}/object-molecules`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ object: testObject.input })
    });
    
    if (!analysisResponse.ok) {
      throw new Error(`Analysis failed: ${analysisResponse.status} ${analysisResponse.statusText}`);
    }
    
    const analysisResult = await analysisResponse.json();
    const molecules = analysisResult.output?.chemicals || [];
    
    console.log(`   ✅ Analysis complete: Found ${molecules.length} molecules`);
    molecules.forEach((mol, i) => {
      console.log(`      ${i + 1}. ${mol.name} (${mol.smiles})`);
    });
    
    if (molecules.length === 0) {
      console.log('   ❌ No molecules found!');
      return false;
    }
    
    // Step 2: Generate SDF files
    console.log('   📁 Generating SDF files...');
    const smilesArray = molecules.map(mol => mol.smiles).filter(Boolean);
    
    const sdfResponse = await fetch(`${baseUrl}/generate-sdfs`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ smiles: smilesArray, overwrite: false })
    });
    
    if (!sdfResponse.ok) {
      throw new Error(`SDF generation failed: ${sdfResponse.status} ${sdfResponse.statusText}`);
    }
    
    const sdfResult = await sdfResponse.json();
    console.log(`   ✅ SDF generation complete: ${sdfResult.sdfPaths?.length || 0} files created`);
    
    if (sdfResult.errors && sdfResult.errors.length > 0) {
      console.log(`   ⚠️  Errors: ${sdfResult.errors.length}`);
      sdfResult.errors.forEach(error => console.log(`      - ${error}`));
    }
    
    // Step 3: Test SDF file accessibility
    console.log('   🔗 Testing SDF file access...');
    let accessibleFiles = 0;
    
    for (const sdfPath of (sdfResult.sdfPaths || [])) {
      try {
        const fileResponse = await fetch(`${baseUrl}${sdfPath}`);
        if (fileResponse.ok) {
          accessibleFiles++;
        }
      } catch (err) {
        console.log(`   ❌ Failed to access ${sdfPath}: ${err.message}`);
      }
    }
    
    console.log(`   ✅ SDF accessibility: ${accessibleFiles}/${sdfResult.sdfPaths?.length || 0} files accessible`);
    
    // Step 4: Create visualization data structure
    const visualizationData = molecules.map((mol, index) => ({
      name: mol.name,
      smiles: mol.smiles,
      sdfData: sdfResult.sdfPaths?.[index] ? `file://${sdfResult.sdfPaths[index]}` : null
    }));
    
    console.log(`   🎯 Visualization data ready for ${visualizationData.length} molecules`);
    
    return true;
    
  } catch (error) {
    console.log(`   ❌ Error: ${error.message}`);
    return false;
  }
}

async function runAllTests() {
  console.log('🚀 Starting Molecular Visualization Tests');
  console.log('==========================================');
  console.log('Make sure the dev server is running (npm run dev)');
  console.log('');
  
  let passed = 0;
  let failed = 0;
  
  for (const testObject of testObjects) {
    const success = await testMolecularAnalysis(testObject);
    if (success) {
      passed++;
    } else {
      failed++;
    }
    
    // Brief pause between tests
    await new Promise(resolve => setTimeout(resolve, 1000));
  }
  
  console.log('\n📊 Test Results');
  console.log('================');
  console.log(`✅ Passed: ${passed}`);
  console.log(`❌ Failed: ${failed}`);
  console.log(`📈 Success Rate: ${((passed / (passed + failed)) * 100).toFixed(1)}%`);
  
  if (failed === 0) {
    console.log('\n🎉 All molecular visualization tests passed!');
    console.log('   Your molecular visualization is working correctly.');
  } else {
    console.log('\n⚠️  Some tests failed. Check the output above for details.');
  }
}

async function runInteractiveTest() {
  const rl = readline.createInterface({
    input: process.stdin,
    output: process.stdout
  });
  
  console.log('🧪 Interactive Molecular Testing');
  console.log('Type an object name to test molecular analysis, or "exit" to quit.');
  console.log('Examples: water, apple, coffee, kale, wine\n');
  
  function askForInput() {
    rl.question('Enter object to analyze: ', async (input) => {
      if (input.toLowerCase() === 'exit') {
        rl.close();
        return;
      }
      
      if (input.trim()) {
        const testObject = {
          name: input,
          input: input,
          expected: ['unknown']
        };
        
        await testMolecularAnalysis(testObject);
      }
      
      console.log('');
      askForInput();
    });
  }
  
  askForInput();
}

// Main execution
const args = process.argv.slice(2);

if (args.includes('--interactive') || args.includes('-i')) {
  runInteractiveTest();
} else if (args.includes('--help') || args.includes('-h')) {
  console.log('Molecular Visualization UI Tester');
  console.log('');
  console.log('Usage:');
  console.log('  node test-molecular-ui.js           # Run all predefined tests');
  console.log('  node test-molecular-ui.js -i        # Interactive testing mode');
  console.log('  node test-molecular-ui.js --help    # Show this help');
  console.log('');
  console.log('Make sure the development server is running before testing.');
} else {
  runAllTests();
}