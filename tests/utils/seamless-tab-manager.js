const puppeteer = require('puppeteer');
const fs = require('fs');
const path = require('path');

class SeamlessTabManager {
  static instance = null;
  
  constructor() {
    if (SeamlessTabManager.instance) {
      return SeamlessTabManager.instance;
    }
    
    this.browser = null;
    this.page = null;
    this.debuggingPort = 9225; // Dedicated port for seamless operation
    this.isConnected = false;
    
    SeamlessTabManager.instance = this;
  }

  static getInstance() {
    if (!SeamlessTabManager.instance) {
      SeamlessTabManager.instance = new SeamlessTabManager();
    }
    return SeamlessTabManager.instance;
  }

  async getPage() {
    if (this.page && this.isConnected) {
      try {
        // Verify page is still responsive
        await this.page.evaluate(() => document.title);
        return this.page;
      } catch (error) {
        // Page is no longer responsive, reconnect
        this.isConnected = false;
      }
    }

    await this.connect();
    return this.page;
  }

  async connect() {
    try {
      // Try connecting to existing browser first
      await this.connectToExisting();
      
      if (!this.browser) {
        // Silently launch new browser if none exists
        await this.launchSilently();
      }

      await this.setupMolecularTab();
      this.isConnected = true;
      
    } catch (error) {
      // Fallback: launch fresh browser
      await this.launchSilently();
      await this.setupMolecularTab();
      this.isConnected = true;
    }
  }

  async connectToExisting() {
    try {
      const browserURL = `http://127.0.0.1:${this.debuggingPort}`;
      this.browser = await puppeteer.connect({ 
        browserURL,
        defaultViewport: { width: 1600, height: 1000 }
      });
      
      // Silent operation - no console logs in normal use
    } catch (error) {
      this.browser = null;
    }
  }

  async launchSilently() {
    this.browser = await puppeteer.launch({
      headless: false, // Visible but silent
      defaultViewport: { width: 1600, height: 1000 },
      args: [
        '--no-sandbox',
        '--disable-setuid-sandbox',
        '--disable-web-security',
        '--no-first-run',
        '--disable-default-apps',
        '--disable-extensions',
        '--disable-background-networking',
        '--disable-background-timer-throttling',
        '--disable-backgrounding-occluded-windows',
        '--disable-renderer-backgrounding',
        '--no-default-browser-check',
        '--disable-infobars',
        `--remote-debugging-port=${this.debuggingPort}`,
        '--user-data-dir=./test/seamless-chrome-profile'
      ],
      slowMo: 100,
      devtools: false // No devtools popup
    });
  }

  async setupMolecularTab() {
    const pages = await this.browser.pages();
    
    // Look for existing molecular tab
    let molecularPage = null;
    for (const page of pages) {
      try {
        const url = page.url();
        if (url.includes('localhost:8080')) {
          molecularPage = page;
          break;
        }
      } catch (error) {
        continue;
      }
    }

    if (molecularPage) {
      this.page = molecularPage;
      // Silently refresh to clean state
      await this.page.goto('http://localhost:8080', { waitUntil: 'networkidle0' });
    } else {
      // Use existing tab or create new one
      this.page = pages.length > 0 ? pages[0] : await this.browser.newPage();
      await this.page.goto('http://localhost:8080', { waitUntil: 'networkidle0' });
    }

    // Silent setup - no console monitoring
    await this.page.waitForSelector('input[type="text"]', { timeout: 10000 });
  }

  async executeTest(testInput) {
    const page = await this.getPage();
    
    // Clear and type input
    await page.evaluate(() => {
      const input = document.querySelector('input[type="text"]');
      if (input) {
        input.value = '';
        input.focus();
        input.dispatchEvent(new Event('input', { bubbles: true }));
      }
    });
    
    await page.type('input[type="text"]', testInput, { delay: 50 });
    await page.keyboard.press('Enter');
    
    return page;
  }

  async takeScreenshot(name) {
    if (!this.page) return null;
    
    const screenshotPath = path.join(__dirname, '../screenshots', `seamless_${name}.png`);
    await this.page.screenshot({ path: screenshotPath, fullPage: true });
    return screenshotPath;
  }

  // Silent cleanup - never closes browser
  cleanup() {
    // Intentionally do nothing - leave browser open for reuse
  }
}

module.exports = SeamlessTabManager;