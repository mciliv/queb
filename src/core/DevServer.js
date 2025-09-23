/**
 * Chrome Auto-Reload Development Server
 * 
 * This module provides automatic Chrome integration for development:
 * - Automatically opens Chrome to the application
 * - Monitors file changes intelligently 
 * - Refreshes only when appropriate changes are made
 * - Minimizes breaking changes and unnecessary reloads
 * 
 * Philosophy: "The best tools are invisible until you need them"
 */

const fs = require('fs');
const path = require('path');
const { spawn } = require('child_process');
const { EventEmitter } = require('events');

class DevServer extends EventEmitter {
  constructor(options = {}) {
    super();
    
    this.config = {
      appUrl: options.appUrl || 'http://localhost:8080',
      chromePort: options.chromePort || 9222,
      debounceMs: options.debounceMs || 500,
      watchDirs: options.watchDirs || [
        'src/client',
        'src/core', 
        'src/server'
      ],
      ignorePatterns: options.ignorePatterns || [
        /node_modules/,
        /\.git/,
        /dist/,
        /logs/,
        /temp/,
        /\.DS_Store/,
        /\.map$/
      ]
    };
    
    this.watchers = new Map();
    this.chromeProcess = null;
    this.isReloading = false;
    this.pendingReloads = new Set();
    this.lastReloadTime = 0;
    this.cdpConnection = null;
    
    // File change tracking
    this.changeQueue = [];
    this.debounceTimer = null;
  }

  /**
   * Start the development server with Chrome integration
   */
  async start() {
    console.log('🚀 Starting Chrome Auto-Reload Development Server...');
    
    try {
      // Start Chrome in debug mode
      await this.startChrome();
      
      // Wait for Chrome to fully load
      await this.waitForChrome();
      
      // Connect to Chrome DevTools Protocol
      await this.connectToCDP();
      
      // Navigate to the application
      await this.navigateToApp();
      
      // Start file watching
      this.startFileWatching();
      
      console.log('✅ Development server ready!');
      console.log(`📱 App: ${this.config.appUrl}`);
      console.log('👀 Watching for file changes...');
      
      this.emit('ready');
      
    } catch (error) {
      console.error('❌ Failed to start development server:', error.message);
      await this.cleanup();
      throw error;
    }
  }

  /**
   * Start Chrome with remote debugging enabled
   */
  async startChrome() {
    // Check if Chrome is already running on the debug port
    if (await this.isChromeRunning()) {
      console.log('🔄 Chrome already running, connecting to existing instance...');
      return;
    }

    const chromeArgs = [
      `--remote-debugging-port=${this.config.chromePort}`,
      '--no-first-run',
      '--no-default-browser-check',
      '--disable-background-timer-throttling',
      '--disable-renderer-backgrounding',
      '--disable-backgrounding-occluded-windows',
      '--disable-ipc-flooding-protection',
      this.config.appUrl
    ];

    // Try different Chrome executable paths
    const chromePaths = this.getChromePaths();
    
    for (const chromePath of chromePaths) {
      try {
        if (fs.existsSync(chromePath)) {
          console.log(`🌐 Starting Chrome: ${chromePath}`);
          this.chromeProcess = spawn(chromePath, chromeArgs, {
            detached: true,
            stdio: 'ignore'
          });
          
          this.chromeProcess.on('error', (error) => {
            console.warn(`⚠️ Chrome process error: ${error.message}`);
          });
          
          return;
        }
      } catch (error) {
        continue; // Try next path
      }
    }
    
    throw new Error('Chrome executable not found. Please install Google Chrome.');
  }

  /**
   * Get possible Chrome executable paths for different OS
   */
  getChromePaths() {
    const platform = process.platform;
    
    if (platform === 'darwin') {
      return [
        '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
        '/Applications/Google Chrome Canary.app/Contents/MacOS/Google Chrome Canary'
      ];
    } else if (platform === 'win32') {
      return [
        'C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe',
        'C:\\Program Files (x86)\\Google\\Chrome\\Application\\chrome.exe',
        `${process.env.LOCALAPPDATA}\\Google\\Chrome\\Application\\chrome.exe`
      ];
    } else {
      return [
        '/usr/bin/google-chrome',
        '/usr/bin/chromium-browser',
        '/usr/bin/chromium',
        '/snap/bin/chromium'
      ];
    }
  }

  /**
   * Check if Chrome is running on the debug port
   */
  async isChromeRunning() {
    try {
      const response = await fetch(`http://localhost:${this.config.chromePort}/json/version`);
      return response.ok;
    } catch (error) {
      return false;
    }
  }

  /**
   * Wait for Chrome to be ready
   */
  async waitForChrome(maxAttempts = 30) {
    for (let i = 0; i < maxAttempts; i++) {
      if (await this.isChromeRunning()) {
        return;
      }
      await new Promise(resolve => setTimeout(resolve, 1000));
    }
    throw new Error('Chrome failed to start within timeout period');
  }

  /**
   * Connect to Chrome DevTools Protocol
   */
  async connectToCDP() {
    try {
      // Get list of tabs
      const response = await fetch(`http://localhost:${this.config.chromePort}/json`);
      const tabs = await response.json();
      
      // Find the tab with our app or create new one
      let targetTab = tabs.find(tab => 
        tab.url && tab.url.includes('localhost') && 
        (tab.url.includes(':8080') || tab.url.includes(':3001'))
      );
      
      if (!targetTab) {
        // Create new tab
        const newTabResponse = await fetch(`http://localhost:${this.config.chromePort}/json/new`);
        targetTab = await newTabResponse.json();
      }
      
      console.log(`🔗 Connected to Chrome tab: ${targetTab.title || 'New Tab'}`);
      this.cdpConnection = targetTab;
      
    } catch (error) {
      throw new Error(`Failed to connect to Chrome: ${error.message}`);
    }
  }

  /**
   * Navigate to the application
   */
  async navigateToApp() {
    if (!this.cdpConnection) return;
    
    try {
      // Use CDP to navigate to our app
      const response = await fetch(`http://localhost:${this.config.chromePort}/json/runtime/evaluate`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          expression: `window.location.href = '${this.config.appUrl}'`
        })
      });
      
      console.log(`🏠 Navigated to: ${this.config.appUrl}`);
      
    } catch (error) {
      console.warn(`⚠️ Navigation warning: ${error.message}`);
    }
  }

  /**
   * Start file watching for intelligent reloads
   */
  startFileWatching() {
    this.config.watchDirs.forEach(dir => {
      const fullPath = path.resolve(process.cwd(), dir);
      
      if (!fs.existsSync(fullPath)) {
        console.warn(`⚠️ Watch directory not found: ${fullPath}`);
        return;
      }
      
      try {
        const watcher = fs.watch(fullPath, { recursive: true }, (eventType, filename) => {
          if (filename) {
            this.handleFileChange(dir, filename, eventType);
          }
        });
        
        this.watchers.set(dir, watcher);
        console.log(`👀 Watching: ${dir}`);
        
      } catch (error) {
        console.warn(`⚠️ Failed to watch ${dir}: ${error.message}`);
      }
    });
  }

  /**
   * Handle file changes intelligently
   */
  handleFileChange(watchDir, filename, eventType) {
    // Ignore patterns
    if (this.config.ignorePatterns.some(pattern => pattern.test(filename))) {
      return;
    }
    
    // Prevent reload loops
    if (this.isReloading) return;
    
    const change = {
      dir: watchDir,
      file: filename,
      type: eventType,
      timestamp: Date.now(),
      category: this.categorizeChange(watchDir, filename)
    };
    
    console.log(`📝 File changed: ${watchDir}/${filename} (${change.category})`);
    
    // Add to change queue
    this.changeQueue.push(change);
    
    // Debounce multiple rapid changes
    if (this.debounceTimer) {
      clearTimeout(this.debounceTimer);
    }
    
    this.debounceTimer = setTimeout(() => {
      this.processChangeQueue();
    }, this.config.debounceMs);
  }

  /**
   * Categorize file changes for appropriate reload strategy
   */
  categorizeChange(watchDir, filename) {
    const ext = path.extname(filename).toLowerCase();
    
    // Frontend changes
    if (watchDir.includes('client') || watchDir.includes('core')) {
      if (['.jsx', '.js', '.css', '.html'].includes(ext)) {
        return 'frontend';
      }
    }
    
    // Backend changes
    if (watchDir.includes('server')) {
      if (['.js', '.json'].includes(ext)) {
        return 'backend';
      }
    }
    
    // Configuration changes
    if (filename.includes('config') || filename.includes('.env')) {
      return 'config';
    }
    
    return 'other';
  }

  /**
   * Process queued changes and determine reload strategy
   */
  async processChangeQueue() {
    if (this.changeQueue.length === 0) return;
    
    const changes = [...this.changeQueue];
    this.changeQueue = [];
    
    // Analyze changes to determine reload strategy
    const categories = new Set(changes.map(c => c.category));
    const hasBackend = categories.has('backend') || categories.has('config');
    const hasFrontend = categories.has('frontend');
    
    try {
      if (hasBackend && hasFrontend) {
        console.log('🔄 Full reload: Backend and frontend changes detected');
        await this.hardReload();
      } else if (hasBackend) {
        console.log('🔄 Page reload: Backend changes detected');
        await this.pageReload();
      } else if (hasFrontend) {
        console.log('🔄 Soft reload: Frontend changes detected');
        await this.softReload();
      } else {
        console.log('ℹ️ No reload needed for these changes');
      }
    } catch (error) {
      console.warn(`⚠️ Reload failed: ${error.message}`);
    }
  }

  /**
   * Soft reload - just refresh the page content
   */
  async softReload() {
    if (!this.cdpConnection) return;
    
    this.isReloading = true;
    
    try {
      await fetch(`http://localhost:${this.config.chromePort}/json/runtime/evaluate`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          expression: 'window.location.reload()'
        })
      });
      
      this.lastReloadTime = Date.now();
      
    } finally {
      setTimeout(() => {
        this.isReloading = false;
      }, 2000);
    }
  }

  /**
   * Page reload - refresh the entire page
   */
  async pageReload() {
    if (!this.cdpConnection) return;
    
    this.isReloading = true;
    
    try {
      await fetch(`http://localhost:${this.config.chromePort}/json/runtime/evaluate`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          expression: 'window.location.reload(true)' // Force reload from server
        })
      });
      
      this.lastReloadTime = Date.now();
      
    } finally {
      setTimeout(() => {
        this.isReloading = false;
      }, 3000);
    }
  }

  /**
   * Hard reload - clear cache and reload
   */
  async hardReload() {
    if (!this.cdpConnection) return;
    
    this.isReloading = true;
    
    try {
      // Clear cache and hard reload
      await fetch(`http://localhost:${this.config.chromePort}/json/runtime/evaluate`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          expression: `
            if ('caches' in window) {
              caches.keys().then(names => names.forEach(name => caches.delete(name)));
            }
            setTimeout(() => window.location.reload(true), 100);
          `
        })
      });
      
      this.lastReloadTime = Date.now();
      
    } finally {
      setTimeout(() => {
        this.isReloading = false;
      }, 5000);
    }
  }

  /**
   * Stop the development server and cleanup
   */
  async stop() {
    console.log('🛑 Stopping development server...');
    await this.cleanup();
    console.log('✅ Development server stopped');
  }

  /**
   * Cleanup resources
   */
  async cleanup() {
    // Stop file watchers
    for (const [dir, watcher] of this.watchers) {
      try {
        watcher.close();
        console.log(`✅ Stopped watching: ${dir}`);
      } catch (error) {
        console.warn(`⚠️ Error stopping watcher for ${dir}: ${error.message}`);
      }
    }
    this.watchers.clear();
    
    // Clear timers
    if (this.debounceTimer) {
      clearTimeout(this.debounceTimer);
    }
    
    // Note: We don't kill Chrome as user might have other tabs open
    // Chrome process will continue running independently
    
    this.emit('stopped');
  }

  /**
   * Get development server status
   */
  getStatus() {
    return {
      running: this.watchers.size > 0,
      chromeConnected: !!this.cdpConnection,
      watchedDirs: Array.from(this.watchers.keys()),
      lastReload: this.lastReloadTime ? new Date(this.lastReloadTime).toISOString() : null,
      appUrl: this.config.appUrl
    };
  }
}

module.exports = DevServer;
