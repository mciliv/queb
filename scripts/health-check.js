#!/usr/bin/env node

const fs = require('fs').promises;
const path = require('path');

async function checkHealth() {
  const healthFile = path.join(__dirname, '../logs/health-status.json');
  const statsFile = path.join(__dirname, '../logs/test-stats.json');
  
  try {
    // Check if health file exists and is recent
    const healthData = await fs.readFile(healthFile, 'utf8');
    const health = JSON.parse(healthData);
    
    const lastCheckTime = new Date(health.lastCheck).getTime();
    const now = Date.now();
    const timeSinceLastCheck = now - lastCheckTime;
    
    // Fail if last check was more than 10 minutes ago
    if (timeSinceLastCheck > 10 * 60 * 1000) {
      console.error('Health check stale - last check was', Math.floor(timeSinceLastCheck / 60000), 'minutes ago');
      process.exit(1);
    }
    
    // Check health status
    if (health.status !== 'healthy') {
      console.error('Health status is', health.status);
      process.exit(1);
    }
    
    // Check test stats
    const statsData = await fs.readFile(statsFile, 'utf8');
    const stats = JSON.parse(statsData);
    
    // Warn if failure rate is high
    if (stats.totalRuns > 10) {
      const failureRate = (stats.failedRuns / stats.totalRuns) * 100;
      if (failureRate > 20) {
        console.error('High test failure rate:', failureRate.toFixed(1) + '%');
        process.exit(1);
      }
    }
    
    console.log('Health check passed ✅');
    console.log('Status:', health.status);
    console.log('Uptime:', Math.floor(health.uptime / 60), 'minutes');
    console.log('Test success rate:', ((stats.successfulRuns / stats.totalRuns) * 100).toFixed(1) + '%');
    
    process.exit(0);
    
  } catch (error) {
    console.error('Health check failed:', error.message);
    process.exit(1);
  }
}

checkHealth();