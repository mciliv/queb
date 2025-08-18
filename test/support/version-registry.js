const fs = require('fs');
const path = require('path');
const { execSync } = require('child_process');

class VersionRegistry {
  constructor() {
    this.registryFile = path.join(__dirname, 'versions.json');
    this.ensureRegistry();
  }

  ensureRegistry() {
    if (!fs.existsSync(this.registryFile)) {
      const defaultRegistry = {
        versions: {},
        currentVersion: 'main',
        lastUpdated: new Date().toISOString()
      };
      fs.writeFileSync(this.registryFile, JSON.stringify(defaultRegistry, null, 2));
    }
  }

  loadRegistry() {
    return JSON.parse(fs.readFileSync(this.registryFile, 'utf8'));
  }

  saveRegistry(registry) {
    registry.lastUpdated = new Date().toISOString();
    fs.writeFileSync(this.registryFile, JSON.stringify(registry, null, 2));
  }

  // Add a new version systematically
  addVersion(versionName, options = {}) {
    const registry = this.loadRegistry();
    
    // Get current git info
    const gitCommit = execSync('git rev-parse HEAD').toString().trim();
    const gitBranch = execSync('git branch --show-current').toString().trim();
    const gitTag = this.getLatestTag();
    
    const versionConfig = {
      name: versionName,
      gitCommit,
      gitBranch,
      gitTag,
      timestamp: new Date().toISOString(),
      description: options.description || `Version ${versionName}`,
      ports: this.allocatePorts(versionName),
      testResults: {},
      ...options
    };

    registry.versions[versionName] = versionConfig;
    this.saveRegistry(registry);
    
    console.log(`✅ Added version: ${versionName}`);
    console.log(`   Commit: ${gitCommit.substring(0, 8)}`);
    console.log(`   Branch: ${gitBranch}`);
    console.log(`   Ports: HTTP=${versionConfig.ports.http}, Debug=${versionConfig.ports.debug}`);
    
    return versionConfig;
  }

  // Systematically allocate ports for a version
  allocatePorts(versionName) {
    const registry = this.loadRegistry();
    const existingVersions = Object.keys(registry.versions);
    
    // Base port calculation: hash version name to get consistent ports
    const hash = this.hashString(versionName);
    const basePort = 3000 + (hash % 100) * 10; // 3000-4000 range
    
    return {
      http: basePort,
      https: basePort + 1,
      debug: 9230 + (hash % 50), // 9230-9280 range
      chrome: 9300 + (hash % 50), // 9300-9350 range
      test: basePort + 2
    };
  }

  hashString(str) {
    let hash = 0;
    for (let i = 0; i < str.length; i++) {
      const char = str.charCodeAt(i);
      hash = ((hash << 5) - hash) + char;
      hash = hash & hash; // Convert to 32-bit integer
    }
    return Math.abs(hash);
  }

  getLatestTag() {
    try {
      return execSync('git describe --tags --abbrev=0').toString().trim();
    } catch (error) {
      return 'no-tags';
    }
  }

  // List all registered versions
  listVersions() {
    const registry = this.loadRegistry();
    return Object.values(registry.versions).map(v => ({
      name: v.name,
      commit: v.gitCommit.substring(0, 8),
      branch: v.gitBranch,
      timestamp: v.timestamp,
      ports: v.ports
    }));
  }

  // Remove a version
  removeVersion(versionName) {
    const registry = this.loadRegistry();
    if (registry.versions[versionName]) {
      delete registry.versions[versionName];
      this.saveRegistry(registry);
      console.log(`✅ Removed version: ${versionName}`);
      return true;
    }
    return false;
  }

  // Get version configuration
  getVersion(versionName) {
    const registry = this.loadRegistry();
    return registry.versions[versionName];
  }

  // Update test results for a version
  updateTestResults(versionName, testResults) {
    const registry = this.loadRegistry();
    if (registry.versions[versionName]) {
      registry.versions[versionName].testResults = {
        ...registry.versions[versionName].testResults,
        ...testResults,
        lastTested: new Date().toISOString()
      };
      this.saveRegistry(registry);
    }
  }

  // Create version from current state
  snapshotCurrent(versionName, description) {
    return this.addVersion(versionName, { 
      description,
      snapshot: true,
      source: 'current'
    });
  }

  // Create version from specific commit
  snapshotFromCommit(versionName, commitHash, description) {
    const currentCommit = execSync('git rev-parse HEAD').toString().trim();
    
    try {
      // Checkout the specific commit temporarily
      execSync(`git checkout ${commitHash}`, { stdio: 'ignore' });
      
      const version = this.addVersion(versionName, {
        description,
        snapshot: true,
        source: 'commit',
        sourceCommit: commitHash
      });
      
      return version;
    } finally {
      // Return to original commit
      execSync(`git checkout ${currentCommit}`, { stdio: 'ignore' });
    }
  }
}

module.exports = VersionRegistry;
