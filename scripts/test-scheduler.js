#!/usr/bin/env node

const cron = require('node-cron');
const fs = require('fs').promises;
const path = require('path');
const { exec } = require('child_process');
const BackgroundTestRunner = require('./background-test-runner');

class TestScheduler {
  constructor() {
    this.runner = new BackgroundTestRunner();
    this.healthCheckInterval = null;
    this.scheduledJobs = [];
    this.statsFile = path.join(__dirname, '../logs/test-stats.json');
    this.healthFile = path.join(__dirname, '../logs/health-status.json');
  }

  log(message, level = 'info') {
    const timestamp = new Date().toISOString();
    const levels = {
      info: '📅',
      success: '✅',
      warning: '⚠️',
      error: '❌',
      health: '💚'
    };
    console.log(`${levels[level]} [${timestamp}] ${message}`);
  }

  async ensureLogDirectory() {
    const logDir = path.join(__dirname, '../logs');
    try {
      await fs.mkdir(logDir, { recursive: true });
    } catch (error) {
      this.log(`Error creating log directory: ${error.message}`, 'error');
    }
  }

  async loadStats() {
    try {
      const data = await fs.readFile(this.statsFile, 'utf8');
      return JSON.parse(data);
    } catch {
      return {
        totalRuns: 0,
        successfulRuns: 0,
        failedRuns: 0,
        lastRun: null,
        averageRunTime: 0
      };
    }
  }

  async saveStats(stats) {
    try {
      await fs.writeFile(this.statsFile, JSON.stringify(stats, null, 2));
    } catch (error) {
      this.log(`Error saving stats: ${error.message}`, 'error');
    }
  }

  async updateHealthStatus(status) {
    const health = {
      status,
      timestamp: new Date().toISOString(),
      uptime: process.uptime(),
      memory: process.memoryUsage(),
      testRunnerActive: this.runner.isRunning,
      lastCheck: new Date().toISOString()
    };

    try {
      await fs.writeFile(this.healthFile, JSON.stringify(health, null, 2));
    } catch (error) {
      this.log(`Error updating health status: ${error.message}`, 'error');
    }
  }

  scheduleJobs() {
    this.log('Initializing test scheduler...', 'info');

    // Every 30 minutes - Quick smoke tests
    const smokeJob = cron.schedule('*/30 * * * *', async () => {
      this.log('Running scheduled smoke tests', 'info');
      const startTime = Date.now();
      
      try {
        await this.runner.runTestSuite('smoke');
        const duration = Date.now() - startTime;
        await this.recordTestRun(true, duration);
      } catch (error) {
        await this.recordTestRun(false, Date.now() - startTime);
        this.log(`Smoke test failed: ${error.message}`, 'error');
      }
    });

    // Every 2 hours - Full backend tests
    const backendJob = cron.schedule('0 */2 * * *', async () => {
      this.log('Running scheduled backend tests', 'info');
      const startTime = Date.now();
      
      try {
        await this.runner.runTestSuite('unit-backend');
        const duration = Date.now() - startTime;
        await this.recordTestRun(true, duration);
      } catch (error) {
        await this.recordTestRun(false, Date.now() - startTime);
        this.log(`Backend test failed: ${error.message}`, 'error');
      }
    });

    // Daily at 3 AM - Full test suite
    const fullTestJob = cron.schedule('0 3 * * *', async () => {
      this.log('Running daily full test suite', 'info');
      const startTime = Date.now();
      
      try {
        await this.runner.runIntelligentTests();
        const duration = Date.now() - startTime;
        await this.recordTestRun(true, duration);
        
        // Clean up old logs after successful run
        await this.cleanupOldLogs();
      } catch (error) {
        await this.recordTestRun(false, Date.now() - startTime);
        this.log(`Full test suite failed: ${error.message}`, 'error');
      }
    });

    // Every 5 minutes - Health check
    const healthJob = cron.schedule('*/5 * * * *', async () => {
      await this.performHealthCheck();
    });

    this.scheduledJobs = [smokeJob, backendJob, fullTestJob, healthJob];
    this.log('Test scheduler initialized with 4 scheduled jobs', 'success');
  }

  async recordTestRun(success, duration) {
    const stats = await this.loadStats();
    stats.totalRuns++;
    if (success) {
      stats.successfulRuns++;
    } else {
      stats.failedRuns++;
    }
    stats.lastRun = new Date().toISOString();
    stats.averageRunTime = ((stats.averageRunTime * (stats.totalRuns - 1)) + duration) / stats.totalRuns;
    
    await this.saveStats(stats);
  }

  async performHealthCheck() {
    const checks = {
      diskSpace: await this.checkDiskSpace(),
      memory: this.checkMemory(),
      gitStatus: await this.checkGitStatus(),
      nodeModules: await this.checkNodeModules()
    };

    const allHealthy = Object.values(checks).every(check => check);
    const status = allHealthy ? 'healthy' : 'degraded';
    
    await this.updateHealthStatus(status);
    
    if (!allHealthy) {
      this.log('Health check detected issues', 'warning');
      Object.entries(checks).forEach(([check, passed]) => {
        if (!passed) {
          this.log(`  ❌ ${check} check failed`, 'warning');
        }
      });
    } else {
      this.log('Health check passed', 'health');
    }
  }

  async checkDiskSpace() {
    return new Promise((resolve) => {
      exec('df -h . | tail -1', (error, stdout) => {
        if (error) {
          resolve(false);
          return;
        }
        const usagePercent = parseInt(stdout.split(/\s+/)[4]);
        resolve(usagePercent < 90); // Fail if > 90% full
      });
    });
  }

  checkMemory() {
    const used = process.memoryUsage();
    const heapPercent = (used.heapUsed / used.heapTotal) * 100;
    return heapPercent < 85; // Fail if > 85% heap used
  }

  async checkGitStatus() {
    return new Promise((resolve) => {
      exec('git status --porcelain', (error, stdout) => {
        if (error) {
          resolve(false);
          return;
        }
        // Warn if too many uncommitted changes
        const changes = stdout.split('\n').filter(line => line.trim());
        resolve(changes.length < 50);
      });
    });
  }

  async checkNodeModules() {
    try {
      await fs.access(path.join(__dirname, '../node_modules'));
      return true;
    } catch {
      return false;
    }
  }

  async cleanupOldLogs() {
    const logDir = path.join(__dirname, '../logs');
    const maxAge = 7 * 24 * 60 * 60 * 1000; // 7 days
    
    try {
      const files = await fs.readdir(logDir);
      const now = Date.now();
      
      for (const file of files) {
        if (file.endsWith('.log')) {
          const filePath = path.join(logDir, file);
          const stats = await fs.stat(filePath);
          
          if (now - stats.mtime.getTime() > maxAge) {
            await fs.unlink(filePath);
            this.log(`Cleaned up old log: ${file}`, 'info');
          }
        }
      }
    } catch (error) {
      this.log(`Error cleaning logs: ${error.message}`, 'error');
    }
  }

  async start() {
    await this.ensureLogDirectory();
    this.scheduleJobs();
    await this.performHealthCheck();
    
    // Initial test run
    this.log('Running initial test suite...', 'info');
    await this.runner.runTestSuite('smoke');
    
    // Graceful shutdown
    process.on('SIGTERM', () => this.shutdown());
    process.on('SIGINT', () => this.shutdown());
    
    this.log('Test scheduler is running. Press Ctrl+C to stop.', 'success');
  }

  shutdown() {
    this.log('Shutting down test scheduler...', 'info');
    this.scheduledJobs.forEach(job => job.stop());
    process.exit(0);
  }
}

// Start scheduler if run directly
if (require.main === module) {
  const scheduler = new TestScheduler();
  scheduler.start().catch(error => {
    console.error('Failed to start scheduler:', error);
    process.exit(1);
  });
}

module.exports = TestScheduler;