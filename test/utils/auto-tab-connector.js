const puppeteer = require('puppeteer');

/**
 * Seamless tab connector that automatically finds or creates 
 * a molecular testing tab without any user intervention
 */
class AutoTabConnector {
  static async getOrCreateTab() {
    const debugPorts = [9225, 9224, 9223, 9222]; // Multiple fallback ports
    
    for (const port of debugPorts) {
      try {
        // Try to connect to existing Chrome instance
        const browser = await puppeteer.connect({ 
          browserURL: `http://127.0.0.1:${port}`,
          defaultViewport: { width: 1600, height: 1000 }
        });
        
        const pages = await browser.pages();
        
        // Look for existing molecular testing tab
        for (const page of pages) {
          try {
            const url = page.url();
            if (url.includes('localhost:8080') || url.includes('localhost:3001')) {
              // Found existing molecular tab, refresh it to clean state
              const targetUrl = process.env.FRONTEND_URL || 'http://localhost:8080';
              await page.goto(targetUrl, { waitUntil: 'networkidle0' });
              await page.waitForSelector('input[type="text"]', { timeout: 10000 });
              
              return { browser, page, reused: true };
            }
          } catch (error) {
            continue;
          }
        }
        
        // Use existing tab and navigate to molecular app
        if (pages.length > 0) {
          const page = pages[0];
          const targetUrl = process.env.FRONTEND_URL || 'http://localhost:8080';
          await page.goto(targetUrl, { waitUntil: 'networkidle0' });
          await page.waitForSelector('input[type="text"]', { timeout: 10000 });
          
          return { browser, page, reused: true };
        }
        
        // Create new tab in existing browser
        const page = await browser.newPage();
        const targetUrl = process.env.FRONTEND_URL || 'http://localhost:8080';
        await page.goto(targetUrl, { waitUntil: 'networkidle0' });
        await page.waitForSelector('input[type="text"]', { timeout: 10000 });
        
        return { browser, page, reused: true };
        
      } catch (error) {
        // Try next port
        continue;
      }
    }
    
    // No existing Chrome found, launch new one silently
    const userDataDir = `/Users/m/mol/test/chrome-molecular-profile-${Date.now()}`;
    const browser = await puppeteer.launch({
      headless: false,
      defaultViewport: { width: 1600, height: 1000 },
      userDataDir: userDataDir,
      args: [
        '--no-sandbox',
        '--disable-setuid-sandbox',
        '--disable-web-security',
        '--no-first-run',
        '--disable-default-apps',
        '--disable-infobars',
        '--remote-debugging-port=9225',
        '--disable-background-timer-throttling',
        '--disable-backgrounding-occluded-windows',
        '--disable-renderer-backgrounding'
      ],
      slowMo: 100,
      devtools: false
    });
    
    const page = await browser.newPage();
    const targetUrl = process.env.FRONTEND_URL || 'http://localhost:8080';
    await page.goto(targetUrl, { waitUntil: 'networkidle0' });
    await page.waitForSelector('input[type="text"]', { timeout: 10000 });
    
    return { browser, page, reused: false };
  }
  
  static async executeTestSilently(testInput) {
    const { page } = await AutoTabConnector.getOrCreateTab();
    
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
}

module.exports = AutoTabConnector;