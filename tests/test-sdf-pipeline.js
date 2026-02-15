const fetch = require('node-fetch');
const fs = require('fs');
const path = require('path');

// Load molecular configuration from centralized config
const molecularConfig = require('../config/molecular-config.json');

async function testSDFPipeline() {
    console.log('🧪 Testing Complete SDF Display Pipeline...\n');

    // Test 1: SDF File Generation
    console.log('1. Testing SDF File Generation');
    // Use configured test molecules instead of hardcoded values
    const testSmiles = Object.values(molecularConfig.TEST_MOLECULES).map(mol => mol.smiles);
    
    try {
        const response = await fetch('http://localhost:8080/generate-sdfs', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ smiles: testSmiles })
        });

        if (response.ok) {
            const result = await response.json();
            console.log(`  ✓ Generated ${result.sdfPaths.length} SDF files`);
            
            // Test 2: SDF File Accessibility
            console.log('\n2. Testing SDF File Accessibility');
            for (const sdfPath of result.sdfPaths) {
                try {
                    const sdfResponse = await fetch(`http://localhost:8080${sdfPath}`);
                    if (sdfResponse.ok) {
                        const sdfContent = await sdfResponse.text();
                        console.log(`  ✓ ${sdfPath} accessible (${sdfContent.length} chars)`);
                        
                        // Validate SDF content
                        if (sdfContent.includes('RDKit') && sdfContent.includes('$$$$')) {
                            console.log(`    ✓ Valid SDF format`);
                        } else {
                            console.log(`    ✗ Invalid SDF format`);
                        }
                    } else {
                        console.log(`  ✗ ${sdfPath} not accessible: ${sdfResponse.status}`);
                    }
                } catch (error) {
                    console.log(`  ✗ ${sdfPath} fetch error: ${error.message}`);
                }
            }
        } else {
            console.log(`  ✗ SDF generation failed: ${response.status}`);
        }
    } catch (error) {
        console.log(`  ✗ API request failed: ${error.message}`);
    }

    // Test 3: File System Check
    console.log('\n3. Testing File System');
    const sdfDir = path.join(__dirname, 'sdf_files');
    if (fs.existsSync(sdfDir)) {
        const files = fs.readdirSync(sdfDir).filter(f => f.endsWith('.sdf'));
        console.log(`  ✓ SDF directory exists with ${files.length} files`);
        
        // Check a few specific files using configured molecules
        const testFiles = testSmiles.map(smiles => `${smiles.replace(/[^a-zA-Z0-9]/g, '_')}.sdf`);
        for (const file of testFiles.slice(0, 2)) { // Test first 2 files
            const filePath = path.join(sdfDir, file);
            if (fs.existsSync(filePath)) {
                const stats = fs.statSync(filePath);
                console.log(`  ✓ ${file} exists (${stats.size} bytes)`);
                
                // Validate SDF content
                const content = fs.readFileSync(filePath, 'utf8');
                if (content.includes('RDKit') && content.includes('$$$$')) {
                    console.log(`    ✓ Valid SDF format`);
                } else {
                    console.log(`    ✗ Invalid SDF format`);
                }
            } else {
                console.log(`  ✗ ${file} missing`);
            }
        }
    } else {
        console.log(`  ✗ SDF directory missing: ${sdfDir}`);
    }

    // Python SDF generator removed; HTTP-based pipeline replaces it

    console.log('\n✅ Pipeline test completed');
}

// Run the test
testSDFPipeline().catch(console.error); 