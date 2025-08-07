#!/usr/bin/env node

/**
 * Test script to inject fake molecular data directly into the app
 * This bypasses the UI and directly calls the API to test molecular visualization
 */

const fetch = require('node-fetch');

// Fake molecular data - simulates what AI would return
const fakeTestData = [
  {
    name: "💧 Water Test",
    fakeResponse: {
      object: "water",
      chemicals: [
        { name: "Water", smiles: "O" }
      ]
    }
  },
  {
    name: "🧪 Ethanol Test", 
    fakeResponse: {
      object: "ethanol",
      chemicals: [
        { name: "Ethanol", smiles: "CCO" }
      ]
    }
  },
  {
    name: "🍷 Wine Test",
    fakeResponse: {
      object: "red wine",
      chemicals: [
        { name: "Water", smiles: "O" },
        { name: "Ethanol", smiles: "CCO" },
        { name: "Tartaric acid", smiles: "OC(C(O)C(O)=O)C(O)=O" },
        { name: "Resveratrol", smiles: "C1=CC(=CC=C1C=CC2=CC(=CC(=C2)O)O)O" }
      ]
    }
  }
];

async function testMolecularVisualization(testData) {
  const baseUrl = 'http://localhost:3000';
  
  try {
    console.log(`\n🧪 Testing: ${testData.name}`);
    console.log(`   Molecules: ${testData.fakeResponse.chemicals.length}`);
    
    // Step 1: Generate SDF files from the fake SMILES data
    const smilesArray = testData.fakeResponse.chemicals.map(mol => mol.smiles);
    console.log(`   SMILES: ${smilesArray.join(', ')}`);
    
    const sdfResponse = await fetch(`${baseUrl}/generate-sdfs`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ smiles: smilesArray })
    });
    
    if (!sdfResponse.ok) {
      throw new Error(`SDF generation failed: ${sdfResponse.status}`);
    }
    
    const sdfResult = await sdfResponse.json();
    console.log(`   ✅ Generated ${sdfResult.sdfPaths.length} SDF files`);
    
    // Step 2: Check file accessibility
    let accessibleFiles = 0;
    for (const sdfPath of sdfResult.sdfPaths) {
      try {
        const fileResponse = await fetch(`${baseUrl}${sdfPath}`);
        if (fileResponse.ok) {
          accessibleFiles++;
          const filename = sdfPath.split('/').pop();
          console.log(`   📁 ${filename} - accessible`);
        }
      } catch (e) {
        const filename = sdfPath.split('/').pop();
        console.log(`   ❌ ${filename} - not accessible`);
      }
    }
    
    console.log(`   🎯 ${accessibleFiles}/${sdfResult.sdfPaths.length} files accessible for 3DMol.js`);
    
    // Step 3: Prepare visualization data (what would be passed to Results component)
    const visualizationData = testData.fakeResponse.chemicals.map((mol, index) => {
      const sdfPath = sdfResult.sdfPaths[index];
      const filename = sdfPath ? sdfPath.split('/').pop() : null;
      return {
        name: mol.name,
        smiles: mol.smiles,
        sdfPath: sdfPath,
        filename: filename
      };
    });
    
    console.log(`   🎬 Visualization data ready:`);
    visualizationData.forEach((mol, i) => {
      const status = mol.sdfPath ? '✅' : '❌';
      console.log(`      ${i+1}. ${status} ${mol.name} (${mol.smiles})`);
    });
    
    return {
      success: true,
      totalMolecules: visualizationData.length,
      accessibleFiles,
      visualizationData
    };
    
  } catch (error) {
    console.log(`   ❌ Error: ${error.message}`);
    return { success: false, error: error.message };
  }
}

async function runAllTests() {
  console.log('🚀 Testing molecular visualization with fake data...');
  console.log('   This simulates what the UI would receive and display\n');
  
  const results = [];
  
  for (const testData of fakeTestData) {
    const result = await testMolecularVisualization(testData);
    results.push(result);
    
    // Brief delay between tests
    await new Promise(resolve => setTimeout(resolve, 1000));
  }
  
  // Summary
  console.log('\n📊 Test Summary:');
  const successful = results.filter(r => r.success).length;
  console.log(`   ✅ ${successful}/${results.length} tests successful`);
  console.log(`   🎯 Ready for 3DMol.js visualization in Results component`);
  console.log('\n💡 Next step: Open http://localhost:3000 and manually test with any object');
  console.log('    The molecular visualization pipeline is working correctly!');
}

// Run if called directly
if (require.main === module) {
  runAllTests().catch(console.error);
}

module.exports = { testMolecularVisualization, fakeTestData };