// test/fixtures/cleanup-registry.js - Central registry for test resource cleanup

class CleanupRegistry {
  constructor() {
    this.resources = new Set();
    this.timers = new Set();
    this.intervals = new Set();
    this.processes = new Set();
  }

  register(resource) {
    this.resources.add(resource);
    return resource;
  }

  registerTimer(timer) {
    this.timers.add(timer);
    return timer;
  }

  registerInterval(interval) {
    this.intervals.add(interval);
    return interval;
  }

  registerProcess(process) {
    this.processes.add(process);
    return process;
  }

  async cleanup() {
    // Clear timers first
    for (const timer of this.timers) {
      clearTimeout(timer);
    }
    this.timers.clear();

    // Clear intervals
    for (const interval of this.intervals) {
      clearInterval(interval);
    }
    this.intervals.clear();

    // Stop processes
    for (const process of this.processes) {
      try {
        if (process.kill) {
          process.kill('SIGTERM');
        }
      } catch (e) {
        console.warn('Process cleanup warning:', e.message);
      }
    }
    this.processes.clear();

    // Cleanup resources
    for (const resource of this.resources) {
      try {
        if (resource && typeof resource.close === 'function') {
          await resource.close();
        } else if (resource && typeof resource.stop === 'function') {
          await resource.stop();
        } else if (resource && typeof resource.destroy === 'function') {
          resource.destroy();
        }
      } catch (e) {
        console.warn('Resource cleanup warning:', e.message);
      }
    }
    this.resources.clear();
  }

  forceCleanup() {
    // Synchronous emergency cleanup
    for (const timer of this.timers) {
      clearTimeout(timer);
    }
    for (const interval of this.intervals) {
      clearInterval(interval);
    }
    for (const process of this.processes) {
      try {
        if (process.kill) process.kill('SIGKILL');
      } catch (e) {
        // Ignore errors in force cleanup
      }
    }
  }
}

// Global registry
global.testCleanupRegistry = global.testCleanupRegistry || new CleanupRegistry();

module.exports = global.testCleanupRegistry;