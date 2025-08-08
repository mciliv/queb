const puppeteer = require('puppeteer');
const fs = require('fs').promises;
const path = require('path');

/**
 * Screenshot Service - Automated Chrome screenshot capture for LLM analysis
 * Creates screenshots that can be read by LLM for visual analysis
 */
class ScreenshotService {
  constructor() {
    this.browser = null;
    this.page = null;
    this.screenshotDir = path.join(__dirname, '../../screenshots');
    this.isInitialized = false;
  }

  async initialize() {
    if (this.isInitialized) return;

    // Ensure screenshot directory exists
    await this.ensureScreenshotDirectory();

    // Launch browser
    this.browser = await puppeteer.launch({
      headless: false, // Keep visible for debugging
      defaultViewport: { width: 1280, height: 720 },
      args: [
        '--no-sandbox',
        '--disable-setuid-sandbox',
        '--disable-web-security',
        '--no-first-run',
        '--disable-default-apps'
      ]
    });

    this.page = await this.browser.newPage();
    this.isInitialized = true;
    console.log('✅ Screenshot service initialized');
  }

  async ensureScreenshotDirectory() {
    try {
      await fs.access(this.screenshotDir);
    } catch {
      await fs.mkdir(this.screenshotDir, { recursive: true });
      console.log(`📁 Created screenshots directory: ${this.screenshotDir}`);
    }
  }

  /**
   * Navigate to the molecular app and take a screenshot
   */
  async captureApp(filename = null) {
    await this.initialize();

    try {
      // Navigate to the app
      await this.page.goto('http://localhost:3000', { 
        waitUntil: 'networkidle0',
        timeout: 15000 
      });

      // Wait for React app to load
      await this.page.waitForSelector('#root', { timeout: 10000 });

      // Generate filename if not provided
      if (!filename) {
        const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
        filename = `app-${timestamp}.png`;
      }

      const screenshotPath = path.join(this.screenshotDir, filename);
      
      // Take screenshot
      await this.page.screenshot({ 
        path: screenshotPath, 
        fullPage: true 
      });

      console.log(`📸 Screenshot saved: ${screenshotPath}`);
      
      return {
        path: screenshotPath,
        filename: filename,
        url: `/api/screenshot/${filename}`,
        timestamp: new Date().toISOString()
      };
    } catch (error) {
      console.error('❌ Screenshot capture failed:', error.message);
      throw error;
    }
  }

  /**
   * Capture screenshot after entering text input
   */
  async captureWithInput(inputText, filename = null) {
    await this.initialize();

    try {
      // Navigate to app
      await this.page.goto('http://localhost:3000', { 
        waitUntil: 'networkidle0',
        timeout: 15000 
      });

      await this.page.waitForSelector('#root', { timeout: 10000 });

      // Find and fill the input
      const textInput = await this.page.$('#object-input');
      if (textInput) {
        await textInput.click();
        await textInput.evaluate(el => el.value = ''); // Clear
        await textInput.type(inputText);
        console.log(`✏️  Entered text: "${inputText}"`);
      }

      // Wait a moment for any UI updates
      await this.page.waitForTimeout(1000);

      // Generate filename if not provided
      if (!filename) {
        const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
        const sanitizedInput = inputText.replace(/[^a-zA-Z0-9]/g, '_').substring(0, 20);
        filename = `input-${sanitizedInput}-${timestamp}.png`;
      }

      const screenshotPath = path.join(this.screenshotDir, filename);
      
      // Take screenshot
      await this.page.screenshot({ 
        path: screenshotPath, 
        fullPage: true 
      });

      console.log(`📸 Screenshot with input saved: ${screenshotPath}`);
      
      return {
        path: screenshotPath,
        filename: filename,
        url: `/api/screenshot/${filename}`,
        inputText: inputText,
        timestamp: new Date().toISOString()
      };
    } catch (error) {
      console.error('❌ Screenshot with input failed:', error.message);
      throw error;
    }
  }

  /**
   * Capture screenshot after triggering analysis
   */
  async captureAnalysis(inputText, filename = null) {
    await this.initialize();

    try {
      // Navigate and input text
      await this.page.goto('http://localhost:3000', { 
        waitUntil: 'networkidle0',
        timeout: 15000 
      });

      await this.page.waitForSelector('#root', { timeout: 10000 });

      // Find input and trigger analysis
      const textInput = await this.page.$('#object-input');
      if (textInput) {
        await textInput.click();
        await textInput.evaluate(el => el.value = '');
        await textInput.type(inputText);
        
        // Press Enter to trigger analysis
        await textInput.press('Enter');
        console.log(`🔬 Triggered analysis for: "${inputText}"`);
        
        // Wait for analysis to complete (adjust timeout as needed)
        await this.page.waitForTimeout(5000);
      }

      // Generate filename
      if (!filename) {
        const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
        const sanitizedInput = inputText.replace(/[^a-zA-Z0-9]/g, '_').substring(0, 20);
        filename = `analysis-${sanitizedInput}-${timestamp}.png`;
      }

      const screenshotPath = path.join(this.screenshotDir, filename);
      
      // Take screenshot
      await this.page.screenshot({ 
        path: screenshotPath, 
        fullPage: true 
      });

      console.log(`📸 Analysis screenshot saved: ${screenshotPath}`);
      
      return {
        path: screenshotPath,
        filename: filename,
        url: `/api/screenshot/${filename}`,
        inputText: inputText,
        timestamp: new Date().toISOString()
      };
    } catch (error) {
      console.error('❌ Analysis screenshot failed:', error.message);
      throw error;
    }
  }

  /**
   * List all available screenshots
   */
  async listScreenshots() {
    await this.ensureScreenshotDirectory();
    
    try {
      const files = await fs.readdir(this.screenshotDir);
      const screenshots = files
        .filter(file => file.endsWith('.png'))
        .map(file => ({
          filename: file,
          path: path.join(this.screenshotDir, file),
          url: `/api/screenshot/${file}`,
          size: null // Could add file stats if needed
        }));

      return screenshots;
    } catch (error) {
      console.error('❌ Failed to list screenshots:', error.message);
      return [];
    }
  }

  /**
   * Get screenshot file path for LLM access
   */
  async getScreenshotPath(filename) {
    const screenshotPath = path.join(this.screenshotDir, filename);
    
    try {
      await fs.access(screenshotPath);
      return screenshotPath;
    } catch {
      throw new Error(`Screenshot not found: ${filename}`);
    }
  }

  /**
   * Clean up old screenshots (keep last N files)
   */
  async cleanupOldScreenshots(keepCount = 10) {
    try {
      const files = await fs.readdir(this.screenshotDir);
      const screenshots = files.filter(file => file.endsWith('.png'));
      
      if (screenshots.length <= keepCount) return;

      // Sort by filename (which includes timestamp)
      screenshots.sort();
      
      // Remove oldest files
      const toDelete = screenshots.slice(0, screenshots.length - keepCount);
      
      for (const file of toDelete) {
        await fs.unlink(path.join(this.screenshotDir, file));
        console.log(`🗑️  Cleaned up old screenshot: ${file}`);
      }
    } catch (error) {
      console.error('❌ Screenshot cleanup failed:', error.message);
    }
  }

  /**
   * Close browser and cleanup
   */
  async close() {
    if (this.browser) {
      await this.browser.close();
      this.browser = null;
      this.page = null;
      this.isInitialized = false;
      console.log('🧹 Screenshot service closed');
    }
  }
}

module.exports = ScreenshotService;
