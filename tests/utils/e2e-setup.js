/**
 * E2E test setup
 * Additional setup for end-to-end tests using Puppeteer
 */

const puppeteer = require('puppeteer');

// Global browser instance
let browser = null;

beforeAll(async () => {
  // Launch browser
  browser = await puppeteer.launch({
    headless: process.env.HEADLESS !== 'false',
    args: [
      '--no-sandbox',
      '--disable-setuid-sandbox',
      '--disable-dev-shm-usage',
      '--disable-accelerated-2d-canvas',
      '--disable-gpu'
    ],
    // Slow down actions for debugging if needed
    slowMo: process.env.SLOW_MO ? parseInt(process.env.SLOW_MO) : 0
  });
  
  global.__BROWSER__ = browser;
  
  console.log('Browser launched for E2E tests');
});

afterAll(async () => {
  // Close browser
  if (browser) {
    await browser.close();
    console.log('Browser closed');
  }
});

// Helper to create a new page with common settings
global.createPage = async (options = {}) => {
  const page = await browser.newPage();
  
  // Set viewport
  await page.setViewport({
    width: options.width || 1280,
    height: options.height || 720
  });
  
  // Set user agent if needed
  if (options.userAgent) {
    await page.setUserAgent(options.userAgent);
  }
  
  // Enable request interception if needed
  if (options.interceptRequests) {
    await page.setRequestInterception(true);
    page.on('request', options.onRequest || ((request) => request.continue()));
  }
  
  // Add console log handler
  page.on('console', (msg) => {
    if (process.env.DEBUG_E2E) {
      console.log(`Browser console [${msg.type()}]:`, msg.text());
    }
  });
  
  // Add error handler
  page.on('pageerror', (error) => {
    console.error('Browser page error:', error.message);
  });
  
  return page;
};

// Helper to wait for element and get it
global.waitForElement = async (page, selector, options = {}) => {
  await page.waitForSelector(selector, {
    visible: options.visible !== false,
    timeout: options.timeout || 10000
  });
  return page.$(selector);
};

// Helper to wait for and click element
global.clickElement = async (page, selector, options = {}) => {
  await global.waitForElement(page, selector, options);
  await page.click(selector, options);
};

// Helper to wait for and type in element
global.typeInElement = async (page, selector, text, options = {}) => {
  await global.waitForElement(page, selector, options);
  
  // Clear existing text if needed
  if (options.clear !== false) {
    await page.click(selector, { clickCount: 3 });
    await page.keyboard.press('Backspace');
  }
  
  await page.type(selector, text, { delay: options.delay || 0 });
};

// Helper to take screenshot on failure
global.takeScreenshotOnFailure = async (page, testName) => {
  if (page && !page.isClosed()) {
    const screenshotPath = `tests/screenshots/failures/${testName}-${Date.now()}.png`;
    try {
      await page.screenshot({ path: screenshotPath, fullPage: true });
      console.log(`Screenshot saved: ${screenshotPath}`);
    } catch (error) {
      console.error('Failed to take screenshot:', error.message);
    }
  }
};

// Helper to wait for network idle
global.waitForNetworkIdle = async (page, options = {}) => {
  await page.waitForLoadState('networkidle', {
    timeout: options.timeout || 30000
  });
};

// Helper to upload file
global.uploadFile = async (page, selector, filePath) => {
  const input = await page.$(selector);
  if (!input) throw new Error(`File input ${selector} not found`);
  await input.uploadFile(filePath);
};

// Add custom matchers for E2E tests
expect.extend({
  async toHaveText(page, selector, expectedText) {
    try {
      await page.waitForSelector(selector, { timeout: 5000 });
      const actualText = await page.$eval(selector, el => el.textContent);
      const pass = actualText.includes(expectedText);
      
      return {
        pass,
        message: () => pass
          ? `Expected element ${selector} not to contain text "${expectedText}"`
          : `Expected element ${selector} to contain text "${expectedText}" but got "${actualText}"`
      };
    } catch (error) {
      return {
        pass: false,
        message: () => `Element ${selector} not found: ${error.message}`
      };
    }
  },
  
  async toBeVisible(page, selector) {
    try {
      await page.waitForSelector(selector, { visible: true, timeout: 5000 });
      return {
        pass: true,
        message: () => `Expected element ${selector} not to be visible`
      };
    } catch (error) {
      return {
        pass: false,
        message: () => `Expected element ${selector} to be visible: ${error.message}`
      };
    }
  }
});
