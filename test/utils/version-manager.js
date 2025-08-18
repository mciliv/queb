const { spawn } = require('child_process');
const puppeteer = require('puppeteer');
const { getVersionConfig, checkVersionAvailability } = require('../support/multi-version-config');

class VersionManager {
  constructor() {
    this.runningVersions = new Map();
    this.browsers = new Map();
  }

  async startVersion(version, options = {}) {
    const config = getVersionConfig(version);
    
    // Check if ports are available
    const availability = await checkVersionAvailability(version);
    if (!availability.available) {
      throw new Error(`Ports not available for version ${version}: ${JSON.stringify(availability.portStatus)}`);
    }

    // Start server process
    const serverProcess = spawn('node', ['backend/api/server.js'], {
      env: {
        ...process.env,
        NODE_ENV: options.nodeEnv || 'development',
        PORT: config.http.toString(),
        HTTPS_PORT: config.https.toString(),
        VERSION_NAME: version,
        NODE_OPTIONS: `--inspect=${config.debug}`
      },
      stdio: options.verbose ? 'inherit' : 'pipe'
    });

    // Wait for server to be ready
    await this.waitForServer(config.http);
    
    this.runningVersions.set(version, {
      config,
      process: serverProcess,
      startTime: Date.now()
    });

    return config;
  }

  async stopVersion(version) {
    const versionInfo = this.runningVersions.get(version);
    if (!versionInfo) {
      throw new Error(`Version ${version} is not running`);
    }

    // Kill server process
    versionInfo.process.kill();
    
    // Close browser if exists
    const browser = this.browsers.get(version);
    if (browser) {
      await browser.close();
      this.browsers.delete(version);
    }

    this.runningVersions.delete(version);
  }

  async connectBrowser(version) {
    const config = getVersionConfig(version);
    
    try {
      // Try connecting to existing Chrome instance
      const browser = await puppeteer.connect({
        browserURL: config.debugUrl,
        defaultViewport: { width: 1600, height: 1000 }
      });
      
      this.browsers.set(version, browser);
      return browser;
      
    } catch (error) {
      // Launch new browser instance
      const browser = await puppeteer.launch({
        headless: false,
        defaultViewport: { width: 1600, height: 1000 },
        userDataDir: `./test/chrome-${version}-profile`,
        args: [
          '--no-sandbox',
          '--disable-setuid-sandbox',
          `--remote-debugging-port=${config.chrome}`,
          '--disable-web-security'
        ]
      });
      
      this.browsers.set(version, browser);
      return browser;
    }
  }

  async runParallelTests(versions, testSuite) {
    const results = new Map();
    
    // Start all versions
    for (const version of versions) {
      try {
        await this.startVersion(version);
        console.log(`✅ ${version} version started on ports: ${JSON.stringify(getVersionConfig(version))}`);
      } catch (error) {
        console.error(`❌ Failed to start ${version}: ${error.message}`);
        results.set(version, { error: error.message });
      }
    }

    // Run tests on each version in parallel
    const testPromises = versions.map(async (version) => {
      if (results.has(version)) return; // Skip if startup failed
      
      try {
        const browser = await this.connectBrowser(version);
        const page = await browser.newPage();
        const config = getVersionConfig(version);
        
        await page.goto(config.baseUrl, { waitUntil: 'networkidle0' });
        
        // Run the test suite
        const testResult = await this.executeTestSuite(page, testSuite);
        results.set(version, { success: true, result: testResult });
        
      } catch (error) {
        results.set(version, { error: error.message });
      }
    });

    await Promise.all(testPromises);
    return results;
  }

  async executeTestSuite(page, testSuite) {
    // Basic molecular test - can be extended
    await page.waitForSelector('input[type="text"]', { timeout: 10000 });
    
    // Test input
    await page.type('input[type="text"]', 'water', { delay: 50 });
    await page.keyboard.press('Enter');
    
    // Wait for results
    await page.waitForSelector('.molecule-container', { timeout: 15000 });
    
    // Capture result
    const molecules = await page.evaluate(() => {
      const containers = document.querySelectorAll('.molecule-container');
      return containers.length;
    });
    
    return { moleculesFound: molecules, testInput: 'water' };
  }

  async waitForServer(port, timeout = 30000) {
    const start = Date.now();
    
    while (Date.now() - start < timeout) {
      try {
        const response = await fetch(`http://localhost:${port}`);
        if (response.ok) return true;
      } catch (error) {
        await new Promise(resolve => setTimeout(resolve, 500));
      }
    }
    
    throw new Error(`Server on port ${port} did not start within ${timeout}ms`);
  }

  getRunningVersions() {
    return Array.from(this.runningVersions.keys());
  }

  async cleanup() {
    for (const version of this.runningVersions.keys()) {
      await this.stopVersion(version);
    }
  }
}

module.exports = VersionManager;
