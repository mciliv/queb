// Simple test runner for component connections
// This can be run from the command line or browser

const { spawn } = require('child_process');
const http = require('http');

class ComponentTestRunner {
  constructor() {
    this.serverUrl = 'http://localhost:8080';
    this.testResults = [];
  }

  // Check if server is running
  async checkServer() {
    return new Promise((resolve) => {
      const req = http.get(this.serverUrl, (res) => {
        resolve(res.statusCode === 200);
      });
      
      req.on('error', () => {
        resolve(false);
      });
      
      req.setTimeout(5000, () => {
        req.destroy();
        resolve(false);
      });
    });
  }

  // Run basic connectivity tests
  async runConnectivityTests() {
    console.log('🔍 Testing server connectivity...');
    
    const serverRunning = await this.checkServer();
    if (!serverRunning) {
      console.log('❌ Server not running. Start with: ./run server');
      return false;
    }
    
    console.log('✅ Server is running');
    return true;
  }

  // Test API endpoints
  async testAPIEndpoints() {
    console.log('\n🌐 Testing API endpoints...');
    
    const endpoints = [
      { path: '/stripe-config', method: 'POST' },
      { path: '/analyze-text', method: 'POST' },
      { path: '/image-molecules', method: 'POST' },
      { path: '/generate-sdfs', method: 'POST' }
    ];
    
    for (const endpoint of endpoints) {
      try {
        const response = await this.makeRequest(endpoint.path, endpoint.method);
        console.log(`✅ ${endpoint.method} ${endpoint.path}: ${response.statusCode}`);
      } catch (error) {
        console.log(`❌ ${endpoint.method} ${endpoint.path}: ${error.message}`);
      }
    }
  }

  // Make HTTP request
  makeRequest(path, method = 'GET', data = null) {
    return new Promise((resolve, reject) => {
      const options = {
        hostname: 'localhost',
        port: 8080,
        path: path,
        method: method,
        headers: {
          'Content-Type': 'application/json'
        }
      };

      const req = http.request(options, (res) => {
        let body = '';
        res.on('data', (chunk) => {
          body += chunk;
        });
        res.on('end', () => {
          resolve({ statusCode: res.statusCode, body });
        });
      });

      req.on('error', reject);
      
      if (data) {
        req.write(JSON.stringify(data));
      }
      
      req.end();
    });
  }

  // Test file structure
  testFileStructure() {
    console.log('\n📁 Testing file structure...');
    
    const fs = require('fs');
    const path = require('path');
    
    const requiredFiles = [
      'frontend/core/index.html',
      'frontend/core/app.js',
      'frontend/components/camera.js',
      'frontend/components/camera-handler.js',
      'frontend/components/ui-utils.js',
      'backend/api/server.js',
      'backend/services/molecular-processor.js'
    ];
    
    let allFilesExist = true;
    
    requiredFiles.forEach(file => {
      const filePath = path.join(process.cwd(), file);
      const exists = fs.existsSync(filePath);
      console.log(`${exists ? '✅' : '❌'} ${file}: ${exists ? 'Found' : 'Missing'}`);
      if (!exists) allFilesExist = false;
    });
    
    return allFilesExist;
  }

  // Test Python dependencies
  async testPythonDependencies() {
    console.log('\n🐍 Testing Python dependencies...');
    
    try {
      const { exec } = require('child_process');
      const util = require('util');
      const execAsync = util.promisify(exec);
      
      // Test Python availability
      const { stdout: pythonVersion } = await execAsync('python --version');
      console.log(`✅ Python: ${pythonVersion.trim()}`);
      
      // Test required Python packages
      const packages = ['rdkit', 'numpy', 'requests'];
      for (const pkg of packages) {
        try {
          await execAsync(`python -c "import ${pkg}"`);
          console.log(`✅ ${pkg}: Available`);
        } catch (error) {
          console.log(`❌ ${pkg}: Missing`);
        }
      }
      
    } catch (error) {
      console.log(`❌ Python test failed: ${error.message}`);
    }
  }

  // Run all tests
  async runAllTests() {
    console.log('🧪 Starting Component Connection Tests\n');
    
    // Test 1: File structure
    const fileStructureOk = this.testFileStructure();
    
    // Test 2: Server connectivity
    const serverOk = await this.runConnectivityTests();
    
    // Test 3: API endpoints
    if (serverOk) {
      await this.testAPIEndpoints();
    }
    
    // Test 4: Python dependencies
    await this.testPythonDependencies();
    
    // Summary
    console.log('\n📊 Test Summary:');
    console.log(`File Structure: ${fileStructureOk ? '✅' : '❌'}`);
    console.log(`Server Connectivity: ${serverOk ? '✅' : '❌'}`);
    
    if (fileStructureOk && serverOk) {
      console.log('\n🎉 Basic component connections are working!');
      console.log('\n📝 Next steps:');
      console.log('1. Open http://localhost:8080 in your browser');
      console.log('2. Open browser console and run: componentTests.runAllTests()');
      console.log('3. Test camera click, text analysis, and photo upload functionality');
    } else {
      console.log('\n⚠️ Some components need attention before testing UI interactions');
    }
  }
}

// Run tests if this file is executed directly
if (require.main === module) {
  const runner = new ComponentTestRunner();
  runner.runAllTests().catch(console.error);
}

module.exports = ComponentTestRunner; 