/**
 * Chrome Development Utilities
 * 
 * Utilities for Chrome DevTools Protocol integration and browser detection
 */

const { spawn } = require('child_process');
const fs = require('fs');

class ChromeUtils {
  /**
   * Check if Chrome is available on the system
   */
  static isAvailable() {
    const chromePaths = ChromeUtils.getChromePaths();
    return chromePaths.some(path => fs.existsSync(path));
  }

  /**
   * Get Chrome executable paths for current platform
   */
  static getChromePaths() {
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
   * Check if Chrome is running with DevTools on specified port
   */
  static async isRunningWithDevTools(port = 9222) {
    try {
      const response = await fetch(`http://localhost:${port}/json/version`);
      if (response.ok) {
        const data = await response.json();
        return {
          running: true,
          version: data.Browser || 'Unknown',
          userAgent: data['User-Agent'] || 'Unknown'
        };
      }
    } catch (error) {
      // Chrome not running or DevTools not enabled
    }
    
    return { running: false };
  }

  /**
   * Get Chrome installation info
   */
  static getInstallationInfo() {
    const chromePaths = ChromeUtils.getChromePaths();
    const available = [];
    
    for (const path of chromePaths) {
      if (fs.existsSync(path)) {
        try {
          const stats = fs.statSync(path);
          available.push({
            path,
            size: stats.size,
            modified: stats.mtime
          });
        } catch (error) {
          // Ignore stat errors
        }
      }
    }
    
    return {
      platform: process.platform,
      available,
      hasChrome: available.length > 0
    };
  }

  /**
   * Launch Chrome with DevTools enabled
   */
  static async launch(options = {}) {
    const {
      port = 9222,
      url = 'about:blank',
      headless = false,
      userDataDir = null
    } = options;

    const chromePaths = ChromeUtils.getChromePaths();
    const availableChrome = chromePaths.find(path => fs.existsSync(path));
    
    if (!availableChrome) {
      throw new Error('Chrome not found. Please install Google Chrome.');
    }

    const args = [
      `--remote-debugging-port=${port}`,
      '--no-first-run',
      '--no-default-browser-check',
      '--disable-background-timer-throttling',
      '--disable-renderer-backgrounding',
      '--disable-backgrounding-occluded-windows',
      '--disable-ipc-flooding-protection'
    ];

    if (headless) {
      args.push('--headless');
    }

    if (userDataDir) {
      args.push(`--user-data-dir=${userDataDir}`);
    }

    if (url !== 'about:blank') {
      args.push(url);
    }

    try {
      const process = spawn(availableChrome, args, {
        detached: true,
        stdio: 'ignore'
      });

      return {
        success: true,
        pid: process.pid,
        port,
        chrome: availableChrome
      };
    } catch (error) {
      throw new Error(`Failed to launch Chrome: ${error.message}`);
    }
  }

  /**
   * Get Chrome DevTools tabs
   */
  static async getTabs(port = 9222) {
    try {
      const response = await fetch(`http://localhost:${port}/json`);
      if (response.ok) {
        return await response.json();
      }
    } catch (error) {
      throw new Error(`Failed to get Chrome tabs: ${error.message}`);
    }
    
    return [];
  }

  /**
   * Create new Chrome tab
   */
  static async createTab(url, port = 9222) {
    try {
      const createUrl = `http://localhost:${port}/json/new`;
      const fullUrl = url ? `${createUrl}?${encodeURIComponent(url)}` : createUrl;
      
      const response = await fetch(fullUrl, { method: 'POST' });
      if (response.ok) {
        return await response.json();
      }
    } catch (error) {
      throw new Error(`Failed to create Chrome tab: ${error.message}`);
    }
    
    throw new Error('Failed to create Chrome tab');
  }

  /**
   * Close Chrome tab
   */
  static async closeTab(tabId, port = 9222) {
    try {
      const response = await fetch(`http://localhost:${port}/json/close/${tabId}`, {
        method: 'POST'
      });
      return response.ok;
    } catch (error) {
      return false;
    }
  }

  /**
   * Execute JavaScript in Chrome tab
   */
  static async executeInTab(tabId, code, port = 9222) {
    try {
      const response = await fetch(`http://localhost:${port}/json/runtime/evaluate`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          expression: code,
          tabId: tabId
        })
      });
      
      if (response.ok) {
        return await response.json();
      }
    } catch (error) {
      throw new Error(`Failed to execute code in Chrome: ${error.message}`);
    }
    
    throw new Error('Failed to execute code in Chrome tab');
  }
}

module.exports = ChromeUtils;
