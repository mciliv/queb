const path = require('path');
const fs = require('fs').promises;

/**
 * Simple error screenshot utility - captures app state when errors occur
 * Minimal, focused functionality without complex API routes
 */
class ErrorScreenshot {
  constructor() {
    this.screenshotDir = path.join(__dirname, '..', 'logs', 'error-screenshots');
    this.isEnabled = process.env.NODE_ENV === 'development';
  }

  /**
   * Capture screenshot when an error occurs
   * Returns screenshot info or null if disabled/failed
   */
  async captureOnError(error, context = {}) {
    if (!this.isEnabled) return null;

    try {
      // Lazy load puppeteer only when needed
      const puppeteer = require('puppeteer');
      
      // Ensure directory exists
      await fs.mkdir(this.screenshotDir, { recursive: true });

      // Launch browser quickly
      const browser = await puppeteer.launch({
        headless: true, // Headless for error capture
        defaultViewport: { width: 1280, height: 720 },
        args: ['--no-sandbox', '--disable-setuid-sandbox']
      });

      const page = await browser.newPage();

      // Navigate to app
      await page.goto('http://localhost:8080', {
        waitUntil: 'domcontentloaded', // Faster than networkidle
        timeout: 5000 // Quick timeout for errors
      });

      // Generate filename with error context
      const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
      const errorType = error.constructor.name || 'Error';
      const filename = `error-${errorType}-${timestamp}.png`;
      const screenshotPath = path.join(this.screenshotDir, filename);

      // Take screenshot
      await page.screenshot({
        path: screenshotPath,
        fullPage: false // Just viewport for speed
      });

      await browser.close();

      console.log(`📸 Error screenshot saved: ${filename}`);
      
      // Cleanup old error screenshots (keep last 5)
      await this.cleanupOldScreenshots(5);

      return {
        path: screenshotPath,
        filename,
        error: error.message,
        context,
        timestamp: new Date().toISOString()
      };
    } catch (screenshotError) {
      // Don't let screenshot errors break the app
      console.warn('⚠️ Error screenshot failed:', screenshotError.message);
      return null;
    }
  }

  /**
   * Simple cleanup of old error screenshots
   */
  async cleanupOldScreenshots(keepCount = 5) {
    try {
      const files = await fs.readdir(this.screenshotDir);
      const screenshots = files
        .filter(file => file.startsWith('error-') && file.endsWith('.png'))
        .sort()
        .reverse(); // Newest first

      if (screenshots.length <= keepCount) return;

      const toDelete = screenshots.slice(keepCount);
      for (const file of toDelete) {
        await fs.unlink(path.join(this.screenshotDir, file));
      }
    } catch (error) {
      // Silent cleanup failure
    }
  }
}

// Export singleton instance
const errorScreenshot = new ErrorScreenshot();

/**
 * Simple function to call when errors occur
 * Usage: await captureErrorScreenshot(error, { userId: 123, action: 'analyze' });
 */
async function captureErrorScreenshot(error, context = {}) {
  return errorScreenshot.captureOnError(error, context);
}

module.exports = {
  captureErrorScreenshot,
  ErrorScreenshot
};
