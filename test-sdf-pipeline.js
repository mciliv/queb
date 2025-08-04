const fetch = require('node-fetch');
const fs = require('fs');
const path = require('path');

async function testSDFPipeline() {
    console.log('🧪 Testing Complete SDF Display Pipeline...\n');

    // Test 1: SDF File Generation
    console.log('1. Testing SDF File Generation');
    const testSmiles = ['CCO', 'O', 'C1=CC=CC=C1'];
    
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
    const sdfDir = path.join(__dirname, 'data', 'sdf_files');
    if (fs.existsSync(sdfDir)) {
        const files = fs.readdirSync(sdfDir).filter(f => f.endsWith('.sdf'));
        console.log(`  ✓ SDF directory exists with ${files.length} files`);
        
        // Check a few specific files
        const testFiles = ['CCO.sdf', 'O.sdf'];
        for (const file of testFiles) {
            const filePath = path.join(sdfDir, file);
            if (fs.existsSync(filePath)) {
                const stats = fs.statSync(filePath);
                console.log(`  ✓ ${file} exists (${stats.size} bytes)`);
            } else {
                console.log(`  ✗ ${file} missing`);
            }
        }
    } else {
        console.log(`  ✗ SDF directory missing: ${sdfDir}`);
    }

    // Test 4: Python SDF Generator
    console.log('\n4. Testing Python SDF Generator');
    const { spawn } = require('child_process');
    
    try {
        const pythonProcess = spawn('python', [
            path.join(__dirname, 'chemistry', 'processors', 'sdf.py'),
            'CCO',
            '--dir',
            'test_sdf_output',
            '--overwrite'
        ], {
            cwd: __dirname
        });

        await new Promise((resolve, reject) => {
            let output = '';
            let error = '';

            pythonProcess.stdout.on('data', (data) => {
                output += data.toString();
            });

            pythonProcess.stderr.on('data', (data) => {
                error += data.toString();
            });

            pythonProcess.on('close', (code) => {
                if (code === 0) {
                    console.log(`  ✓ Python SDF generator successful`);
                    if (output) console.log(`    Output: ${output.trim()}`);
                } else {
                    console.log(`  ✗ Python SDF generator failed (code ${code})`);
                    if (error) console.log(`    Error: ${error.trim()}`);
                }
                resolve();
            });
        });
    } catch (error) {
        console.log(`  ✗ Python process error: ${error.message}`);
    }

    console.log('\n✅ Pipeline test completed');
}

// Run the test
testSDFPipeline().catch(console.error); 