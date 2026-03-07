const { spawn, exec } = require('child_process');
const fs = require('fs');
const path = require('path');

class AppPlacement {
  constructor(options = {}) {
    this.projectRoot = options.projectRoot || process.cwd();
    
    // Load configuration from file
    const configPath = path.join(this.projectRoot, 'config/app-placement.js');
    let baseConfig = {};
    
    try {
      baseConfig = require(configPath);
    } catch (error) {
      console.log('⚠️ Using default app placement config (config file not found)');
      baseConfig = this.getDefaultConfig();
    }
    
    // Merge with runtime options
    this.config = {
      ...baseConfig,
      // URLs
      httpUrl: `http://localhost:${options.httpPort || baseConfig.urls?.desktop?.defaultPort || 8080}`,
      httpsUrl: `https://localhost:${options.httpsPort || baseConfig.urls?.mobile?.defaultPort || 3001}`,
      
      // Update paths to be absolute
      singleTab: {
        ...baseConfig.singleTab,
        userDataDir: path.join(this.projectRoot, baseConfig.singleTab?.userDataDir || '.chrome-profile')
      }
    };
  }

  getDefaultConfig() {
    return {
      singleTab: { enabled: true, debugPort: 9222, userDataDir: '.chrome-profile' },
      urls: { mobile: { defaultPort: 3001 }, desktop: { defaultPort: 8080 } },
      window: { desktop: { width: 1400, height: 900 } }
    };
  }

  /**
   * Get the optimal URL based on device and requirements
   * HTTP-first strategy for better dev experience
   */
  getOptimalUrl(forMobile = false) {
    const isMobileDevice = forMobile || this.isMobile();
    
    // Desktop: Always use HTTP for better dev experience
    if (!isMobileDevice || this.config.urls?.desktop?.forceHttp) {
      console.log('🚀 Using HTTP for better dev experience');
      return this.config.httpUrl;
    }
    
    // Mobile: Only use HTTPS if camera access required
    if (isMobileDevice && this.config.urls?.mobile?.requiresCamera) {
      console.log('📱 Using HTTPS for mobile camera access');
      return this.config.httpsUrl;
    }
    
    // Default to HTTP for best dev experience
    console.log('🔧 Defaulting to HTTP for optimal development');
    return this.config.httpUrl;
  }

  /**
   * Get Chrome launch arguments for single-tab strategy
   */
  getChromeArgs(targetUrl) {
    const chromeConfig = this.config.singleTab?.chromeArgs || {};
    const args = [];

    // Add configured argument groups
    if (chromeConfig.security) args.push(...chromeConfig.security);
    if (chromeConfig.performance) args.push(...chromeConfig.performance);
    if (chromeConfig.ui) args.push(...chromeConfig.ui);

    // Add debugging and user data dir
    args.push(`--remote-debugging-port=${this.config.singleTab?.debugPort || 9222}`);
    args.push(`--user-data-dir=${this.config.singleTab?.userDataDir}`);

    // Add window size for desktop
    if (!this.isMobile()) {
      const windowConfig = this.config.window?.desktop || { width: 1400, height: 900 };
      args.push(`--window-size=${windowConfig.width},${windowConfig.height}`);
      
      if (windowConfig.devTools) {
        args.push('--auto-open-devtools-for-tabs');
      }
    }
    
    return args;
  }

  /**
   * Detect if running on mobile device
   */
  isMobile() {
    const detection = this.config.deviceDetection || {};
    
    // Check for forced override
    if (detection.forceMobile) return true;
    if (detection.forceDesktop) return false;
    
    // Check user agent
    const userAgent = process.env.USER_AGENT || '';
    const mobilePatterns = detection.mobileUserAgents || ['Mobile', 'Android', 'iPhone', 'iPad'];
    
    return mobilePatterns.some(pattern => userAgent.includes(pattern));
  }

  /**
   * Clean up existing Chrome processes and data
   */
  async cleanup() {
    try {
      // Kill Chrome processes with our debug port
      exec(`pkill -f "chrome.*--remote-debugging-port=${this.config.singleTab.debugPort}" || true`, { stdio: 'ignore' });
      
      // Clean user data directory if it exists
      if (fs.existsSync(this.config.singleTab.userDataDir)) {
        const { rimraf } = require('rimraf');
        await rimraf(this.config.singleTab.userDataDir);
      }
      
      console.log('✅ Cleaned up Chrome processes and data');
    } catch (error) {
      console.log(`⚠️ Cleanup warning: ${error.message}`);
    }
  }

  /**
   * Launch or connect to single Chrome tab
   */
  async launchSingleTab(targetUrl = null) {
    const url = targetUrl || this.getOptimalUrl(this.isMobile());
    
    try {
      // Try to connect to existing Chrome first
      const existing = await this.connectToExistingTab(url);
      if (existing) {
        console.log('🔄 Reusing existing Chrome tab');
        return existing;
      }
    } catch (error) {
      console.log('🆕 No existing Chrome found, launching new instance');
    }

    // Launch new Chrome instance
    return this.launchNewChrome(url);
  }

  /**
   * Connect to existing Chrome tab
   */
  async connectToExistingTab(targetUrl) {
    const puppeteer = require('puppeteer');
    
    try {
      const browser = await puppeteer.connect({
        browserURL: `http://127.0.0.1:${this.config.singleTab.debugPort}`,
        defaultViewport: null
      });

      const pages = await browser.pages();
      
      // Look for existing molecular tab
      for (const page of pages) {
        const url = page.url();
        if (url.includes('localhost:8080') || url.includes('localhost:3001')) {
          // Found existing tab, navigate to target URL
          await page.goto(targetUrl, { waitUntil: 'networkidle0' });
          return { browser, page, reused: true };
        }
      }

      // Use first available page
      if (pages.length > 0) {
        const page = pages[0];
        await page.goto(targetUrl, { waitUntil: 'networkidle0' });
        return { browser, page, reused: true };
      }

      return null;
    } catch (error) {
      return null;
    }
  }

  /**
   * Launch new Chrome instance with single-tab strategy
   */
  async launchNewChrome(targetUrl) {
    const args = this.getChromeArgs(targetUrl);
    
    // Ensure user data directory exists
    if (!fs.existsSync(this.config.singleTab.userDataDir)) {
      fs.mkdirSync(this.config.singleTab.userDataDir, { recursive: true });
    }

    const platform = process.platform;
    let chromePath;

    // Determine Chrome executable path
    if (platform === 'darwin') {
      chromePath = '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome';
    } else if (platform === 'win32') {
      chromePath = 'chrome';
    } else {
      chromePath = 'google-chrome';
    }

    // Launch Chrome
    const chromeProcess = spawn(chromePath, [...args, targetUrl], {
      detached: true,
      stdio: 'ignore'
    });

    chromeProcess.unref();

    // Wait a moment for Chrome to start
    await new Promise(resolve => setTimeout(resolve, 2000));

    // Connect to the launched Chrome
    return this.connectToExistingTab(targetUrl);
  }

  /**
   * Open app with fallback strategy
   */
  async openApp(options = {}) {
    const targetUrl = options.url || this.getOptimalUrl();
    
    console.log(`🌐 Opening app: ${targetUrl}`);
    console.log(`📱 Mobile optimized: ${this.isMobile() ? 'Yes' : 'No'}`);
    console.log(`🔐 HTTPS required: ${targetUrl.startsWith('https') ? 'Yes' : 'No'}`);

    try {
      if (this.config.singleTab.enabled) {
        const result = await this.launchSingleTab(targetUrl);
        if (result) {
          console.log('✅ Single tab strategy successful');
          return result;
        }
      }
    } catch (error) {
      console.log(`⚠️ Single tab strategy failed: ${error.message}`);
    }

    // Fallback to system default browser
    console.log('🔄 Falling back to system browser');
    this.openSystemBrowser(targetUrl);
    return null;
  }

  /**
   * Open in system default browser (fallback)
   */
  openSystemBrowser(url) {
    const platform = process.platform;
    
    try {
      if (platform === 'darwin') {
        exec(`open "${url}"`, { stdio: 'ignore' });
      } else if (platform === 'win32') {
        exec(`start "${url}"`, { shell: true, stdio: 'ignore' });
      } else {
        exec(`xdg-open "${url}"`, { stdio: 'ignore' });
      }
      console.log('🌐 Opened in system browser');
    } catch (error) {
      console.log(`❌ Could not open browser: ${error.message}`);
      console.log(`💡 Manual: ${url}`);
    }
  }

  /**
   * Get placement configuration
   */
  getConfig() {
    return { ...this.config };
  }

  /**
   * Update placement configuration
   */
  updateConfig(updates) {
    this.config = { ...this.config, ...updates };
  }
}

module.exports = AppPlacement;
