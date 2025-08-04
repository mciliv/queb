// connection-health.js - Frontend component for displaying connection health status

class ConnectionHealthMonitor {
  constructor() {
    this.isMonitoring = false;
    this.checkInterval = null;
    this.lastHealthData = null;
  }

  // Initialize connection health monitoring
  initialize() {
    this.createHealthDisplay();
    this.bindEvents();
    
    // Auto-check health on app load (dev mode only)
    if (location.hostname === 'localhost' || location.hostname === '127.0.0.1') {
      setTimeout(() => this.checkConnectionHealth(), 2000);
    }
  }

  // Create health display UI
  createHealthDisplay() {
    // Only show in development mode
    if (location.hostname !== 'localhost' && location.hostname !== '127.0.0.1') {
      return;
    }

    const healthContainer = document.createElement('div');
    healthContainer.id = 'connection-health';
    healthContainer.className = 'connection-health hidden';
    healthContainer.innerHTML = `
      <div class="health-header">
        <span class="health-title">🔗 Connection Health</span>
        <button id="health-refresh" class="health-btn">↻</button>
        <button id="health-toggle" class="health-btn">✕</button>
      </div>
      <div class="health-content">
        <div id="health-summary" class="health-summary">
          <span class="health-percentage">--</span>
          <span class="health-status">Unknown</span>
        </div>
        <div id="health-details" class="health-details"></div>
      </div>
    `;

    // Add CSS
    const style = document.createElement('style');
    style.textContent = `
      .connection-health {
        position: fixed;
        top: 10px;
        right: 10px;
        width: 300px;
        background: rgba(0, 0, 0, 0.9);
        color: white;
        border: 1px solid #333;
        border-radius: 4px;
        font-family: monospace;
        font-size: 11px;
        z-index: 10000;
        transition: all 0.3s ease;
      }
      
      .connection-health.hidden {
        transform: translateX(320px);
      }
      
      .health-header {
        display: flex;
        justify-content: space-between;
        align-items: center;
        padding: 8px;
        background: rgba(255, 255, 255, 0.1);
        border-bottom: 1px solid #333;
      }
      
      .health-title {
        font-weight: bold;
      }
      
      .health-btn {
        background: none;
        border: none;
        color: white;
        cursor: pointer;
        padding: 2px 6px;
        margin-left: 4px;
        border-radius: 2px;
      }
      
      .health-btn:hover {
        background: rgba(255, 255, 255, 0.2);
      }
      
      .health-content {
        padding: 8px;
      }
      
      .health-summary {
        display: flex;
        justify-content: space-between;
        margin-bottom: 8px;
        font-weight: bold;
      }
      
      .health-percentage {
        font-size: 14px;
      }
      
      .health-status.excellent { color: #4CAF50; }
      .health-status.good { color: #FFC107; }
      .health-status.critical { color: #F44336; }
      
      .health-details {
        font-size: 10px;
        line-height: 1.3;
      }
      
      .health-step {
        display: flex;
        justify-content: space-between;
        margin: 2px 0;
      }
      
      .health-step.pass { color: #4CAF50; }
      .health-step.partial { color: #FFC107; }
      .health-step.fail { color: #F44336; }
      
      .health-loading {
        text-align: center;
        color: #888;
      }
    `;
    
    document.head.appendChild(style);
    document.body.appendChild(healthContainer);
  }

  // Bind event handlers
  bindEvents() {
    const refreshBtn = document.getElementById('health-refresh');
    const toggleBtn = document.getElementById('health-toggle');
    
    if (refreshBtn) {
      refreshBtn.addEventListener('click', () => this.checkConnectionHealth());
    }
    
    if (toggleBtn) {
      toggleBtn.addEventListener('click', () => this.toggleHealthDisplay());
    }

    // Keyboard shortcut: Ctrl+Shift+H to toggle health display
    document.addEventListener('keydown', (e) => {
      if (e.ctrlKey && e.shiftKey && e.key === 'H') {
        e.preventDefault();
        this.toggleHealthDisplay();
      }
    });
  }

  // Toggle health display visibility
  toggleHealthDisplay() {
    const healthContainer = document.getElementById('connection-health');
    if (healthContainer) {
      healthContainer.classList.toggle('hidden');
      
      // If showing and no data, check health
      if (!healthContainer.classList.contains('hidden') && !this.lastHealthData) {
        this.checkConnectionHealth();
      }
    }
  }

  // Check connection health from backend
  async checkConnectionHealth() {
    const healthSummary = document.getElementById('health-summary');
    const healthDetails = document.getElementById('health-details');
    
    if (!healthSummary || !healthDetails) return;
    
    // Show loading state
    healthDetails.innerHTML = '<div class="health-loading">Checking connections...</div>';
    
    try {
      const response = await fetch('/health/connections');
      
      if (!response.ok) {
        throw new Error(`Health check failed: ${response.status}`);
      }
      
      const healthData = await response.json();
      this.lastHealthData = healthData;
      this.displayHealthData(healthData);
      
    } catch (error) {
      console.error('Failed to check connection health:', error);
      this.displayHealthError(error.message);
    }
  }

  // Display health data in UI
  displayHealthData(healthData) {
    const healthSummary = document.getElementById('health-summary');
    const healthDetails = document.getElementById('health-details');
    
    if (!healthSummary || !healthDetails) return;
    
    // Update summary
    const percentage = Math.round(healthData.connectionHealth);
    const statusClass = healthData.summary || 'unknown';
    
    healthSummary.innerHTML = `
      <span class="health-percentage">${percentage}%</span>
      <span class="health-status ${statusClass}">${statusClass.toUpperCase()}</span>
    `;
    
    // Update details
    const steps = [
      { name: 'Frontend', key: 'step5_frontend' },
      { name: 'SDF Files', key: 'step4_sdfServing' },
      { name: 'SDF Gen', key: 'step3_sdfGeneration' },
      { name: 'Analysis', key: 'step2_molecularAnalysis' },
      { name: 'Input', key: 'step1_inputValidation' }
    ];
    
    const detailsHTML = steps.map(step => {
      const result = healthData.details[step.key];
      const icon = result.status === 'pass' ? '✅' : 
                   result.status === 'partial' ? '⚠️' : '❌';
      
      return `
        <div class="health-step ${result.status}">
          <span>${icon} ${step.name}</span>
          <span>${result.status.toUpperCase()}</span>
        </div>
      `;
    }).join('');
    
    healthDetails.innerHTML = detailsHTML;
  }

  // Display health check error
  displayHealthError(errorMessage) {
    const healthSummary = document.getElementById('health-summary');
    const healthDetails = document.getElementById('health-details');
    
    if (!healthSummary || !healthDetails) return;
    
    healthSummary.innerHTML = `
      <span class="health-percentage">ERR</span>
      <span class="health-status critical">ERROR</span>
    `;
    
    healthDetails.innerHTML = `
      <div class="health-step fail">
        <span>❌ Health Check Failed</span>
      </div>
      <div style="margin-top: 4px; font-size: 9px; color: #888;">
        ${errorMessage}
      </div>
    `;
  }

  // Auto-refresh health data
  startAutoRefresh(intervalMs = 30000) {
    if (this.checkInterval) {
      clearInterval(this.checkInterval);
    }
    
    this.checkInterval = setInterval(() => {
      const healthContainer = document.getElementById('connection-health');
      if (healthContainer && !healthContainer.classList.contains('hidden')) {
        this.checkConnectionHealth();
      }
    }, intervalMs);
  }

  // Stop auto-refresh
  stopAutoRefresh() {
    if (this.checkInterval) {
      clearInterval(this.checkInterval);
      this.checkInterval = null;
    }
  }
}

// Global instance
window.connectionHealthMonitor = new ConnectionHealthMonitor();

// Auto-initialize in development mode
if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', () => {
    window.connectionHealthMonitor.initialize();
  });
} else {
  window.connectionHealthMonitor.initialize();
}

export { ConnectionHealthMonitor };