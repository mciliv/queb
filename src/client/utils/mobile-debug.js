// Mobile Debug Console
// Lightweight debugging tool for mobile devices

class MobileDebug {
  constructor() {
    this.isEnabled = false;
    this.logs = [];
    this.maxLogs = 100;
    this.panel = null;
    this.isVisible = false;
  }

  /**
   * Initialize mobile debug console
   */
  init() {
    // Only enable in development or if debug param is present
    const isDev = process.env.NODE_ENV === 'development';
    const hasDebugParam = new URLSearchParams(window.location.search).has('debug');
    const isLocalhost = window.location.hostname === 'localhost' || window.location.hostname === '127.0.0.1';
    
    if (isDev || hasDebugParam || isLocalhost) {
      this.enable();
    }
  }

  /**
   * Enable mobile debugging
   */
  enable() {
    if (this.isEnabled) return;
    
    this.isEnabled = true;
    this.createPanel();
    this.interceptConsole();
    this.addTouchGestures();
    
    // Add debug info to page
    this.log('🔧 Mobile Debug Console Active');
    this.log(`📱 User Agent: ${navigator.userAgent}`);
    this.log(`📏 Screen: ${screen.width}x${screen.height}`);
    this.log(`🌐 Viewport: ${window.innerWidth}x${window.innerHeight}`);
    this.log(`🔋 Connection: ${navigator.connection?.effectiveType || 'unknown'}`);
    
    console.log('📱 Mobile Debug Console initialized. Triple-tap to toggle.');
  }

  /**
   * Create debug panel
   */
  createPanel() {
    this.panel = document.createElement('div');
    this.panel.id = 'mobile-debug-panel';
    this.panel.innerHTML = `
      <div class="debug-header">
        <span>🔧 Debug Console</span>
        <button class="debug-close" onclick="window.mobileDebug.toggle()">×</button>
      </div>
      <div class="debug-content">
        <div class="debug-logs"></div>
        <div class="debug-controls">
          <button onclick="window.mobileDebug.clear()">Clear</button>
          <button onclick="window.mobileDebug.exportLogs()">Export</button>
          <button onclick="window.mobileDebug.testCamera()">Test Camera</button>
        </div>
      </div>
    `;

    // Add styles
    const style = document.createElement('style');
    style.textContent = `
      #mobile-debug-panel {
        position: fixed;
        top: 20px;
        right: 20px;
        width: 300px;
        max-height: 400px;
        background: rgba(0, 0, 0, 0.95);
        color: #00ff00;
        font-family: 'Courier New', monospace;
        font-size: 12px;
        border: 1px solid #333;
        border-radius: 8px;
        z-index: 10000;
        display: none;
        backdrop-filter: blur(10px);
      }
      
      .debug-header {
        padding: 8px 12px;
        background: rgba(255, 255, 255, 0.1);
        display: flex;
        justify-content: space-between;
        align-items: center;
        font-weight: bold;
      }
      
      .debug-close {
        background: none;
        border: none;
        color: #ff6666;
        font-size: 16px;
        cursor: pointer;
        padding: 0;
        width: 20px;
        height: 20px;
      }
      
      .debug-content {
        padding: 12px;
      }
      
      .debug-logs {
        max-height: 250px;
        overflow-y: auto;
        margin-bottom: 12px;
        padding: 8px;
        background: rgba(255, 255, 255, 0.05);
        border-radius: 4px;
      }
      
      .debug-log {
        margin-bottom: 4px;
        word-wrap: break-word;
      }
      
      .debug-log.error { color: #ff6666; }
      .debug-log.warn { color: #ffaa00; }
      .debug-log.info { color: #66aaff; }
      
      .debug-controls {
        display: flex;
        gap: 8px;
      }
      
      .debug-controls button {
        flex: 1;
        padding: 6px;
        background: rgba(255, 255, 255, 0.1);
        border: 1px solid #333;
        color: #fff;
        border-radius: 4px;
        font-size: 10px;
        cursor: pointer;
      }
      
      .debug-controls button:active {
        background: rgba(255, 255, 255, 0.2);
      }
      
      @media (max-width: 480px) {
        #mobile-debug-panel {
          width: calc(100vw - 40px);
          right: 20px;
        }
      }
    `;
    
    document.head.appendChild(style);
    document.body.appendChild(this.panel);
  }

  /**
   * Intercept console methods
   */
  interceptConsole() {
    const originalLog = console.log;
    const originalError = console.error;
    const originalWarn = console.warn;
    const originalInfo = console.info;

    console.log = (...args) => {
      originalLog.apply(console, args);
      this.log(args.join(' '), 'log');
    };

    console.error = (...args) => {
      originalError.apply(console, args);
      this.log(args.join(' '), 'error');
    };

    console.warn = (...args) => {
      originalWarn.apply(console, args);
      this.log(args.join(' '), 'warn');
    };

    console.info = (...args) => {
      originalInfo.apply(console, args);
      this.log(args.join(' '), 'info');
    };
  }

  /**
   * Add touch gestures for toggle
   */
  addTouchGestures() {
    let tapCount = 0;
    let tapTimer = null;

    document.addEventListener('touchend', (e) => {
      tapCount++;
      
      if (tapCount === 1) {
        tapTimer = setTimeout(() => {
          tapCount = 0;
        }, 500);
      } else if (tapCount === 3) {
        clearTimeout(tapTimer);
        tapCount = 0;
        this.toggle();
      }
    });
  }

  /**
   * Log message
   */
  log(message, type = 'log') {
    const timestamp = new Date().toLocaleTimeString();
    const logEntry = {
      message,
      type,
      timestamp
    };
    
    this.logs.push(logEntry);
    
    if (this.logs.length > this.maxLogs) {
      this.logs.shift();
    }
    
    this.updatePanel();
  }

  /**
   * Update debug panel
   */
  updatePanel() {
    if (!this.panel || !this.isVisible) return;
    
    const logsContainer = this.panel.querySelector('.debug-logs');
    logsContainer.innerHTML = this.logs
      .map(log => `<div class="debug-log ${log.type}">[${log.timestamp}] ${log.message}</div>`)
      .join('');
    
    logsContainer.scrollTop = logsContainer.scrollHeight;
  }

  /**
   * Toggle debug panel visibility
   */
  toggle() {
    this.isVisible = !this.isVisible;
    this.panel.style.display = this.isVisible ? 'block' : 'none';
    
    if (this.isVisible) {
      this.updatePanel();
    }
  }

  /**
   * Clear logs
   */
  clear() {
    this.logs = [];
    this.updatePanel();
  }

  /**
   * Export logs
   */
  exportLogs() {
    const logText = this.logs
      .map(log => `[${log.timestamp}] ${log.type.toUpperCase()}: ${log.message}`)
      .join('\n');
    
    const blob = new Blob([logText], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `debug-logs-${new Date().toISOString().slice(0, 19)}.txt`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
    
    this.log('📄 Logs exported');
  }

  /**
   * Test camera functionality
   */
  async testCamera() {
    this.log('📷 Testing camera...');
    
    try {
      const stream = await navigator.mediaDevices.getUserMedia({ video: true });
      this.log('✅ Camera access granted');
      
      const videoTrack = stream.getVideoTracks()[0];
      const capabilities = videoTrack.getCapabilities();
      const settings = videoTrack.getSettings();
      
      this.log(`📹 Camera: ${settings.width}x${settings.height}`);
      this.log(`🎯 Facing: ${settings.facingMode || 'unknown'}`);
      this.log(`⚡ FPS: ${settings.frameRate || 'unknown'}`);
      
      if (capabilities.focusMode) {
        this.log(`🔍 Focus modes: ${capabilities.focusMode.join(', ')}`);
      }
      
      stream.getTracks().forEach(track => track.stop());
      this.log('✅ Camera test complete');
      
    } catch (error) {
      this.log(`❌ Camera error: ${error.message}`, 'error');
    }
  }
}

// Initialize global instance
window.mobileDebug = new MobileDebug();

// Auto-initialize
if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', () => window.mobileDebug.init());
} else {
  window.mobileDebug.init();
}

export default MobileDebug;

