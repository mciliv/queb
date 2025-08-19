const fetch = require('node-fetch');

async function testAPIEndpoint() {
    console.log('🧪 Testing /generate-sdfs API endpoint...\n');

    const testCases = [
        { smiles: ['CCO'], name: 'Single SMILES' },
        { smiles: ['CCO', 'O', 'C1=CC=CC=C1'], name: 'Multiple SMILES' },
        { smiles: ['INVALID_SMILES'], name: 'Invalid SMILES' }
    ];

    for (const testCase of testCases) {
        console.log(`Testing: ${testCase.name}`);
        console.log(`SMILES: ${testCase.smiles.join(', ')}`);
        
        try {
            const response = await fetch('http://localhost:8080/generate-sdfs', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ smiles: testCase.smiles })
            });

            console.log(`  Status: ${response.status} ${response.statusText}`);

            if (response.ok) {
                const result = await response.json();
                console.log(`  Response:`, JSON.stringify(result, null, 2));
                
                if (result.sdfPaths) {
                    console.log(`  ✓ SDF paths: ${result.sdfPaths.length}`);
                    for (const path of result.sdfPaths) {
                        console.log(`    - ${path}`);
                    }
                }
                
                if (result.errors && result.errors.length > 0) {
                    console.log(`  ⚠ Errors: ${result.errors.length}`);
                    for (const error of result.errors) {
                        console.log(`    - ${error}`);
                    }
                }
                
                if (result.skipped && result.skipped.length > 0) {
                    console.log(`  ⏭ Skipped: ${result.skipped.length}`);
                    for (const skipped of result.skipped) {
                        console.log(`    - ${skipped}`);
                    }
                }
            } else {
                const errorText = await response.text();
                console.log(`  ✗ Error response: ${errorText}`);
            }
        } catch (error) {
            console.log(`  ✗ Request failed: ${error.message}`);
        }
        
        console.log('');
    }
}

// Run the test
testAPIEndpoint().catch(console.error); 