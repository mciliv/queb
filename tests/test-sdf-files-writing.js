#!/usr/bin/env node
/**
 * SDF Files Writing Test
 * 
 * This script tests that SDF files are properly written to the tests/sdf_files directory
 * and validates the file system operations, content format, and HTTP accessibility.
 */

const fs = require('fs');
const path = require('path');
const { spawn } = require('child_process');

// Test configuration
const TEST_SDF_DIR = path.join(__dirname, 'sdf_files');
const TEST_MOLECULES = [
  { name: 'ethanol', smiles: 'CCO', expectedFile: 'CCO.sdf' },
  { name: 'water', smiles: 'O', expectedFile: 'O.sdf' },
  { name: 'caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', expectedFile: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C.sdf' },
  { name: 'benzene', smiles: 'C1=CC=CC=C1', expectedFile: 'C1=CC=CC=C1.sdf' }
];

// Ensure test directory exists
if (!fs.existsSync(TEST_SDF_DIR)) {
  fs.mkdirSync(TEST_SDF_DIR, { recursive: true });
  console.log(`✅ Created test directory: ${TEST_SDF_DIR}`);
}

/**
 * Validate SDF file format
 */
function validateSDFContent(filePath) {
  try {
    const content = fs.readFileSync(filePath, 'utf8');
    
    // Basic SDF format checks
    const hasEndMarker = content.includes('M  END') || content.includes('$$$$');
    const hasReasonableLength = content.length > 50;
    const hasAtomBlock = /^\s*\d+\s+\d+\s+\d+/m.test(content);
    
    return {
      valid: hasEndMarker && hasReasonableLength,
      hasEndMarker,
      hasReasonableLength,
      hasAtomBlock,
      fileSize: fs.statSync(filePath).size,
      content: content.substring(0, 200) + (content.length > 200 ? '...' : '')
    };
  } catch (error) {
    return {
      valid: false,
      error: error.message
    };
  }
}

/**
 * Test SDF file generation via API
 */
async function testSDFGeneration() {
  console.log('🧪 Testing SDF File Generation via API...\n');
  
  const testSmiles = TEST_MOLECULES.map(mol => mol.smiles);
  
  try {
    const response = await fetch('http://localhost:8080/generate-sdfs', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ smiles: testSmiles })
    });

    if (response.ok) {
      const result = await response.json();
      console.log(`✅ Generated ${result.sdfPaths.length} SDF files`);
      console.log(`   Errors: ${result.errors.length}`);
      console.log(`   Skipped: ${result.skipped.length}`);
      
      return result.sdfPaths;
    } else {
      console.log(`❌ SDF generation failed: ${response.status}`);
      return [];
    }
  } catch (error) {
    console.log(`❌ API request failed: ${error.message}`);
    return [];
  }
}

/**
 * Test file system operations
 */
function testFileSystemOperations() {
  console.log('\n📁 Testing File System Operations...\n');
  
  if (!fs.existsSync(TEST_SDF_DIR)) {
    console.log(`❌ SDF directory missing: ${TEST_SDF_DIR}`);
    return false;
  }
  
  const files = fs.readdirSync(TEST_SDF_DIR).filter(f => f.endsWith('.sdf'));
  console.log(`✅ SDF directory exists with ${files.length} files`);
  
  // Test each expected file
  let validFiles = 0;
  for (const molecule of TEST_MOLECULES) {
    const filePath = path.join(TEST_SDF_DIR, molecule.expectedFile);
    
    if (fs.existsSync(filePath)) {
      const validation = validateSDFContent(filePath);
      if (validation.valid) {
        console.log(`✅ ${molecule.name} (${molecule.smiles}) - Valid SDF (${validation.fileSize} bytes)`);
        validFiles++;
      } else {
        console.log(`❌ ${molecule.name} (${molecule.smiles}) - Invalid SDF format`);
        console.log(`   Issues: ${JSON.stringify(validation)}`);
      }
    } else {
      console.log(`⚠️  ${molecule.name} (${molecule.smiles}) - File not found: ${molecule.expectedFile}`);
    }
  }
  
  console.log(`\n📊 File System Summary: ${validFiles}/${TEST_MOLECULES.length} valid files`);
  return validFiles > 0;
}

/**
 * Test HTTP accessibility
 */
async function testHTTPAccessibility(sdfPaths) {
  console.log('\n🌐 Testing HTTP Accessibility...\n');
  
  let accessibleFiles = 0;
  
  for (const sdfPath of sdfPaths) {
    try {
      const response = await fetch(`http://localhost:8080${sdfPath}`);
      
      if (response.ok) {
        const content = await response.text();
        const validation = validateSDFContent(null, content);
        
        if (validation.valid) {
          console.log(`✅ ${sdfPath} - Accessible and valid (${content.length} chars)`);
          accessibleFiles++;
        } else {
          console.log(`❌ ${sdfPath} - Accessible but invalid format`);
        }
      } else {
        console.log(`❌ ${sdfPath} - Not accessible: ${response.status}`);
      }
    } catch (error) {
      console.log(`❌ ${sdfPath} - Fetch error: ${error.message}`);
    }
  }
  
  console.log(`\n📊 HTTP Summary: ${accessibleFiles}/${sdfPaths.length} accessible files`);
  return accessibleFiles > 0;
}

/**
 * Test Python SDF script
 */
function testPythonSDFScript() {
  console.log('\n🐍 Testing Python SDF Script...\n');
  
  return new Promise((resolve) => {
    const pythonScript = path.join(__dirname, '../src/server/molecular-docking-research/sdf.py');
    const args = [pythonScript, 'CCO', '--dir', TEST_SDF_DIR];
    
    const child = spawn('python', args, {
      stdio: ['pipe', 'pipe', 'pipe'],
      timeout: 30000
    });
    
    let stdout = '';
    let stderr = '';
    
    child.stdout.on('data', (data) => {
      stdout += data.toString();
    });
    
    child.stderr.on('data', (data) => {
      stderr += data.toString();
    });
    
    child.on('close', (code) => {
      if (code === 0) {
        console.log('✅ Python SDF script executed successfully');
        console.log(`   Output: ${stdout.trim()}`);
        
        // Check if file was created
        const expectedFile = path.join(TEST_SDF_DIR, 'CCO.sdf');
        if (fs.existsSync(expectedFile)) {
          const validation = validateSDFContent(expectedFile);
          if (validation.valid) {
            console.log(`✅ Generated file is valid SDF (${validation.fileSize} bytes)`);
            resolve(true);
          } else {
            console.log(`❌ Generated file is invalid SDF format`);
            resolve(false);
          }
        } else {
          console.log(`❌ Expected file not created: ${expectedFile}`);
          resolve(false);
        }
      } else {
        console.log(`❌ Python SDF script failed with code ${code}`);
        console.log(`   Error: ${stderr.trim()}`);
        resolve(false);
      }
    });
    
    child.on('error', (error) => {
      console.log(`❌ Python SDF script error: ${error.message}`);
      resolve(false);
    });
  });
}

/**
 * Main test function
 */
async function runTests() {
  console.log('🧪 SDF Files Writing Test Suite');
  console.log('================================\n');
  
  console.log(`📁 Test Directory: ${TEST_SDF_DIR}`);
  console.log(`🧬 Test Molecules: ${TEST_MOLECULES.length}\n`);
  
  const results = {
    apiGeneration: false,
    fileSystem: false,
    httpAccess: false,
    pythonScript: false
  };
  
  // Test 1: API Generation
  const sdfPaths = await testSDFGeneration();
  results.apiGeneration = sdfPaths.length > 0;
  
  // Test 2: File System Operations
  results.fileSystem = testFileSystemOperations();
  
  // Test 3: HTTP Accessibility
  if (sdfPaths.length > 0) {
    results.httpAccess = await testHTTPAccessibility(sdfPaths);
  }
  
  // Test 4: Python Script
  results.pythonScript = await testPythonSDFScript();
  
  // Summary
  console.log('\n📊 Test Results Summary');
  console.log('======================');
  console.log(`API Generation: ${results.apiGeneration ? '✅ PASS' : '❌ FAIL'}`);
  console.log(`File System: ${results.fileSystem ? '✅ PASS' : '❌ FAIL'}`);
  console.log(`HTTP Access: ${results.httpAccess ? '✅ PASS' : '❌ FAIL'}`);
  console.log(`Python Script: ${results.pythonScript ? '✅ PASS' : '❌ FAIL'}`);
  
  const passedTests = Object.values(results).filter(Boolean).length;
  const totalTests = Object.keys(results).length;
  
  console.log(`\n🎯 Overall: ${passedTests}/${totalTests} tests passed`);
  
  if (passedTests === totalTests) {
    console.log('🎉 All tests passed! SDF files are being written correctly to tests/sdf_files');
  } else {
    console.log('⚠️  Some tests failed. Check the output above for details.');
  }
  
  return passedTests === totalTests;
}

// Run tests if this script is executed directly
if (require.main === module) {
  runTests().then(success => {
    process.exit(success ? 0 : 1);
  }).catch(error => {
    console.error('💥 Test suite failed:', error);
    process.exit(1);
  });
}

module.exports = { runTests, validateSDFContent, TEST_SDF_DIR };