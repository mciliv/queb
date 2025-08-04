const fs = require('fs');
const path = require('path');
const { spawn } = require('child_process');

async function testSDFGeneration() {
    console.log('🧪 Testing SDF file generation...\n');

    const testCases = [
        { smiles: 'CCO', name: 'Ethanol' },
        { smiles: 'O', name: 'Water' },
        { smiles: 'C1=CC=CC=C1', name: 'Benzene' },
        { smiles: 'CC(=O)O', name: 'Acetic Acid' }
    ];

    const testDir = 'test_sdf_output';
    if (!fs.existsSync(testDir)) {
        fs.mkdirSync(testDir);
    }

    for (const testCase of testCases) {
        console.log(`Testing: ${testCase.name} (${testCase.smiles})`);
        
        try {
            // Generate SDF file
            const result = await generateSDF(testCase.smiles, testDir);
            
            if (result.success) {
                console.log(`  ✓ SDF generated: ${result.filepath}`);
                
                // Validate SDF content
                const validation = validateSDFContent(result.filepath, testCase.smiles);
                if (validation.valid) {
                    console.log(`  ✓ SDF content valid`);
                    console.log(`    - Atoms: ${validation.atomCount}`);
                    console.log(`    - Bonds: ${validation.bondCount}`);
                    console.log(`    - SMILES property: ${validation.hasSmilesProperty}`);
                } else {
                    console.log(`  ✗ SDF content invalid: ${validation.error}`);
                }
            } else {
                console.log(`  ✗ SDF generation failed: ${result.error}`);
            }
        } catch (error) {
            console.log(`  ✗ Test failed: ${error.message}`);
        }
        
        console.log('');
    }

    // Test existing SDF files
    console.log('🧪 Testing existing SDF files...\n');
    const existingSDFs = [
        'data/sdf_files/CCO.sdf',
        'test/fixtures/CCO.sdf'
    ];

    for (const sdfPath of existingSDFs) {
        if (fs.existsSync(sdfPath)) {
            console.log(`Testing existing file: ${sdfPath}`);
            const validation = validateSDFContent(sdfPath);
            if (validation.valid) {
                console.log(`  ✓ File valid`);
                console.log(`    - Atoms: ${validation.atomCount}`);
                console.log(`    - Bonds: ${validation.bondCount}`);
                console.log(`    - SMILES property: ${validation.hasSmilesProperty}`);
            } else {
                console.log(`  ✗ File invalid: ${validation.error}`);
            }
        } else {
            console.log(`  ✗ File not found: ${sdfPath}`);
        }
        console.log('');
    }
}

function generateSDF(smiles, outputDir) {
    return new Promise((resolve) => {
        const pythonProcess = spawn('python', [
            path.join(__dirname, 'chemistry', 'processors', 'sdf.py'),
            smiles,
            '--dir',
            outputDir,
            '--overwrite'
        ], {
            cwd: __dirname
        });

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
                const filepath = path.join(outputDir, `${smiles}.sdf`);
                if (fs.existsSync(filepath)) {
                    resolve({ success: true, filepath });
                } else {
                    resolve({ success: false, error: 'File not created' });
                }
            } else {
                resolve({ success: false, error: `Process failed with code ${code}: ${error}` });
            }
        });
    });
}

function validateSDFContent(filepath, expectedSmiles = null) {
    try {
        const content = fs.readFileSync(filepath, 'utf8');
        const lines = content.split('\n');
        
        // Check for basic SDF structure
        if (lines.length < 4) {
            return { valid: false, error: 'File too short' };
        }

        // Check header line (should contain RDKit) - skip empty lines
        let headerLine = 0;
        while (headerLine < lines.length && lines[headerLine].trim() === '') {
            headerLine++;
        }
        if (headerLine >= lines.length || !lines[headerLine].includes('RDKit')) {
            return { valid: false, error: 'Missing RDKit header' };
        }

        // Parse counts line (after header)
        const countsLine = lines[headerLine + 2];
        const counts = countsLine.trim().split(/\s+/);
        if (counts.length < 3) {
            return { valid: false, error: 'Invalid counts line' };
        }

        const atomCount = parseInt(counts[0]);
        const bondCount = parseInt(counts[1]);

        if (isNaN(atomCount) || isNaN(bondCount)) {
            return { valid: false, error: 'Invalid atom/bond counts' };
        }

        // Check for M END
        const hasMend = lines.some(line => line.trim() === 'M  END');
        if (!hasMend) {
            return { valid: false, error: 'Missing M END' };
        }

        // Check for SMILES property
        const hasSmilesProperty = content.includes('>  <SMILES>') || content.includes('> <SMILES>');

        // Check for $$$$ terminator
        const hasTerminator = content.includes('$$$$');

        if (!hasTerminator) {
            return { valid: false, error: 'Missing $$$$ terminator' };
        }

        return {
            valid: true,
            atomCount,
            bondCount,
            hasSmilesProperty,
            hasTerminator
        };

    } catch (error) {
        return { valid: false, error: error.message };
    }
}

// Run the test
testSDFGeneration().catch(console.error); 