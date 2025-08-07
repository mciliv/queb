const puppeteer = require('puppeteer');
const fs = require('fs');
const path = require('path');

class ChromeTabManager {
  constructor() {
    this.browser = null;
    this.page = null;
    this.tabSessionFile = path.join(__dirname, '../.tab-session.json');
    this.debuggingPort = 9224; // Dedicated port for molecular testing
  }

  async getOrCreateTab() {
    try {
      // Try to connect to existing Chrome instance
      await this.connectToExistingBrowser();
      
      if (!this.browser) {
        // Launch new browser if none exists
        await this.launchNewBrowser();
      }

      // Get or create the molecular testing tab
      await this.getOrCreateMolecularTab();
      
      return this.page;
    } catch (error) {
      console.log(`⚠️  Tab manager error: ${error.message}`);
      // Fallback: launch fresh browser
      await this.launchNewBrowser();
      await this.getOrCreateMolecularTab();
      return this.page;
    }
  }

  async connectToExistingBrowser() {
    try {
      const browserURL = `http://127.0.0.1:${this.debuggingPort}`;
      this.browser = await puppeteer.connect({ 
        browserURL,
        defaultViewport: { width: 1600, height: 1000 }
      });
      console.log('✅ Connected to existing Chrome instance');
      
      // Verify browser is still responsive
      const pages = await this.browser.pages();
      console.log(`   📄 Found ${pages.length} existing tabs`);
      
    } catch (error) {
      console.log(`   ℹ️  No existing Chrome instance found: ${error.message}`);
      this.browser = null;
    }
  }

  async launchNewBrowser() {
    console.log('🚀 Launching new Chrome instance...');
    
    this.browser = await puppeteer.launch({
      headless: false,
      defaultViewport: { width: 1600, height: 1000 },
      args: [
        '--no-sandbox',
        '--disable-setuid-sandbox',
        '--disable-web-security',
        '--start-maximized',
        '--no-first-run',
        '--disable-default-apps',
        '--disable-extensions',
        `--remote-debugging-port=${this.debuggingPort}`,
        '--user-data-dir=./test/chrome-molecular-profile' // Persistent profile
      ],
      slowMo: 200,
      ignoreDefaultArgs: ['--disable-extensions'] // Allow extensions if needed
    });

    console.log(`✅ Chrome launched with debugging on port ${this.debuggingPort}`);
    
    // Save session info
    await this.saveSessionInfo();
  }

  async getOrCreateMolecularTab() {
    const pages = await this.browser.pages();
    
    // Look for existing molecular testing tab
    let molecularPage = null;
    
    for (const page of pages) {
      try {
        const url = page.url();
        if (url.includes('localhost:3001') || url.includes('molecular')) {
          molecularPage = page;
          console.log('♻️  Found existing molecular testing tab');
          break;
        }
      } catch (error) {
        // Skip pages that can't be accessed
        continue;
      }
    }

    if (molecularPage) {
      this.page = molecularPage;
      
      // Refresh the page to ensure clean state
      console.log('🔄 Refreshing molecular tab...');
      await this.page.goto('http://localhost:3001', { waitUntil: 'networkidle0' });
      
    } else {
      // Create new tab for molecular testing
      console.log('📄 Creating new molecular testing tab...');
      
      if (pages.length > 0) {
        // Use existing blank tab if available
        this.page = pages[0];
      } else {
        // Create completely new tab
        this.page = await this.browser.newPage();
      }
      
      // Navigate to molecular app
      await this.page.goto('http://localhost:3001', { waitUntil: 'networkidle0' });
    }

    // Set up page monitoring
    await this.setupPageMonitoring();
    
    // Wait for app to be ready
    await this.waitForMolecularApp();
    
    console.log('✅ Molecular testing tab ready!');
  }

  async setupPageMonitoring() {
    // Enhanced console logging for molecular analysis
    this.page.on('console', msg => {
      const text = msg.text();
      if (text.includes('🧪') || text.includes('🔗') || text.includes('📨') || 
          text.includes('API') || text.includes('molecular') || text.includes('SMILES')) {
        console.log(`   📱 App: ${text}`);
      }
    });

    // Monitor for page crashes or navigation away
    this.page.on('error', error => {
      console.log(`⚠️  Page error: ${error.message}`);
    });

    this.page.on('pageerror', error => {
      console.log(`⚠️  Page JavaScript error: ${error.message}`);
    });
  }

  async waitForMolecularApp() {
    try {
      // Wait for the molecular app to be ready
      await this.page.waitForSelector('input[type="text"]', { timeout: 15000 });
      await new Promise(resolve => setTimeout(resolve, 2000));
      
      // Verify app is interactive
      const inputExists = await this.page.$('input[type="text"]');
      if (!inputExists) {
        throw new Error('Molecular app input not found');
      }
      
    } catch (error) {
      console.log(`⚠️  Molecular app not ready: ${error.message}`);
      throw error;
    }
  }

  async saveSessionInfo() {
    const sessionInfo = {
      debuggingPort: this.debuggingPort,
      timestamp: new Date().toISOString(),
      pid: process.pid
    };
    
    try {
      fs.writeFileSync(this.tabSessionFile, JSON.stringify(sessionInfo, null, 2));
    } catch (error) {
      console.log(`⚠️  Could not save session info: ${error.message}`);
    }
  }

  async clearInput() {
    if (!this.page) return;
    
    await this.page.evaluate(() => {
      const input = document.querySelector('input[type="text"]');
      if (input) {
        input.value = '';
        input.focus();
        input.dispatchEvent(new Event('input', { bubbles: true }));
      }
    });
    
    await new Promise(resolve => setTimeout(resolve, 300));
  }

  async typeInput(text) {
    if (!this.page) return;
    
    // Clear and type in one operation to avoid timing issues
    await this.page.evaluate((inputText) => {
      const input = document.querySelector('input[type="text"]');
      if (input) {
        input.value = '';
        input.focus();
        input.value = inputText;
        input.dispatchEvent(new Event('input', { bubbles: true }));
        input.dispatchEvent(new Event('change', { bubbles: true }));
      }
    }, text);
    
    await new Promise(resolve => setTimeout(resolve, 200));
  }

  async triggerAnalysis() {
    if (!this.page) return;
    
    await this.page.keyboard.press('Enter');
  }

  async takeScreenshot(name) {
    if (!this.page) return null;
    
    const screenshotPath = path.join(__dirname, '../screenshots', `persistent_tab_${name}.png`);
    await this.page.screenshot({ path: screenshotPath, fullPage: true });
    return screenshotPath;
  }

  async cleanup() {
    // Don't close browser - leave it open for inspection
    console.log('💡 Leaving Chrome open for manual inspection');
    console.log(`   Debug port: ${this.debuggingPort}`);
    console.log('   Close browser manually when done');
  }

  async getTabInfo() {
    if (!this.page) return null;
    
    const url = this.page.url();
    const title = await this.page.title();
    
    return { url, title };
  }
}

module.exports = ChromeTabManager;