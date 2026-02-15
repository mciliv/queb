const http = require('http');

async function testSDFHTTP() {
    console.log('🌐 Testing SDF file HTTP access...\n');

    const testUrls = [
        'http://localhost:8080/sdf_files/CCO.sdf',
        'http://localhost:8080/test/support/fixtures/CCO.sdf'
    ];

    for (const url of testUrls) {
        try {
            console.log(`Testing: ${url}`);
            
            const response = await fetch(url);
            
            if (response.ok) {
                const content = await response.text();
                console.log(`  ✓ HTTP 200 OK`);
                console.log(`  - Content length: ${content.length} characters`);
                console.log(`  - Has RDKit header: ${content.includes('RDKit')}`);
                console.log(`  - Has SMILES property: ${content.includes('>  <SMILES>') || content.includes('> <SMILES>')}`);
                console.log(`  - Has $$$$ terminator: ${content.includes('$$$$')}`);
                
                // Check first few lines
                const lines = content.split('\n');
                console.log(`  - First 3 non-empty lines:`);
                let lineCount = 0;
                for (let i = 0; i < lines.length && lineCount < 3; i++) {
                    if (lines[i].trim() !== '') {
                        console.log(`    ${lineCount + 1}: "${lines[i].trim()}"`);
                        lineCount++;
                    }
                }
            } else {
                console.log(`  ✗ HTTP ${response.status}: ${response.statusText}`);
            }
        } catch (error) {
            console.log(`  ✗ Error: ${error.message}`);
        }
        
        console.log('');
    }
}

// Wait a bit for server to start, then test
setTimeout(testSDFHTTP, 2000); 