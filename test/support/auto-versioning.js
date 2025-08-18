const { execSync } = require('child_process');
const VersionRegistry = require('./version-registry');

class AutoVersioning {
  constructor() {
    this.registry = new VersionRegistry();
  }

  // Parse commit message for version type
  parseCommitMessage(message) {
    const msg = message.toLowerCase();
    
    // Semantic versioning patterns
    if (msg.includes('breaking') || msg.includes('major') || msg.includes('!:')) {
      return 'major';
    }
    if (msg.includes('feat') || msg.includes('feature') || msg.includes('add')) {
      return 'minor';
    }
    if (msg.includes('fix') || msg.includes('patch') || msg.includes('bug')) {
      return 'patch';
    }
    
    // Default to patch for any change
    return 'patch';
  }

  // Get commit-based version name
  getCommitVersion(commitHash) {
    const shortHash = commitHash.substring(0, 8);
    const commitNumber = this.getCommitNumber(commitHash);
    return `c${commitNumber}-${shortHash}`;
  }

  // Get sequential commit number
  getCommitNumber(commitHash) {
    try {
      const result = execSync(`git rev-list --count ${commitHash}`).toString().trim();
      return parseInt(result);
    } catch (error) {
      return 0;
    }
  }

  // Get semantic version from git history
  getSemanticVersion() {
    try {
      // Try to get latest tag
      const latestTag = execSync('git describe --tags --abbrev=0 2>/dev/null || echo "v0.0.0"').toString().trim();
      
      // Count commits since last tag
      const commitsSinceTag = execSync(`git rev-list --count ${latestTag}..HEAD 2>/dev/null || echo "0"`).toString().trim();
      
      if (commitsSinceTag === '0') {
        return latestTag;
      }
      
      // Parse version number
      const versionMatch = latestTag.match(/v?(\d+)\.(\d+)\.(\d+)/);
      if (!versionMatch) {
        return `v0.0.${commitsSinceTag}`;
      }
      
      let [, major, minor, patch] = versionMatch.map(Number);
      
      // Get commits since tag and determine version bump
      const commits = this.getCommitsSince(latestTag);
      const versionType = this.analyzeCommits(commits);
      
      switch (versionType) {
        case 'major':
          major++;
          minor = 0;
          patch = 0;
          break;
        case 'minor':
          minor++;
          patch = 0;
          break;
        case 'patch':
        default:
          patch++;
          break;
      }
      
      return `v${major}.${minor}.${patch}`;
      
    } catch (error) {
      return 'v0.0.1';
    }
  }

  getCommitsSince(tag) {
    try {
      const result = execSync(`git log ${tag}..HEAD --oneline`).toString().trim();
      return result ? result.split('\n') : [];
    } catch (error) {
      return [];
    }
  }

  analyzeCommits(commits) {
    let hasBreaking = false;
    let hasFeature = false;
    
    for (const commit of commits) {
      const type = this.parseCommitMessage(commit);
      if (type === 'major') hasBreaking = true;
      if (type === 'minor') hasFeature = true;
    }
    
    if (hasBreaking) return 'major';
    if (hasFeature) return 'minor';
    return 'patch';
  }

  // Auto-add current version
  addCurrentVersion(options = {}) {
    const currentCommit = execSync('git rev-parse HEAD').toString().trim();
    const commitMessage = execSync('git log -1 --pretty=%B').toString().trim();
    
    let versionName;
    if (options.useCommitNumber) {
      versionName = this.getCommitVersion(currentCommit);
    } else if (options.useSemantic) {
      versionName = this.getSemanticVersion();
    } else {
      // Default: use both
      const semantic = this.getSemanticVersion();
      const commit = this.getCommitVersion(currentCommit);
      versionName = `${semantic}-${commit}`;
    }
    
    return this.registry.addVersion(versionName, {
      description: commitMessage,
      auto: true,
      versionType: this.parseCommitMessage(commitMessage),
      ...options
    });
  }

  // Auto-add versions from git history
  addVersionsFromHistory(count = 5, options = {}) {
    const commits = execSync(`git log --oneline -${count}`).toString().trim().split('\n');
    const versions = [];
    
    commits.forEach((commit, index) => {
      const [hash, ...messageParts] = commit.split(' ');
      const message = messageParts.join(' ');
      
      let versionName;
      if (options.useCommitNumber) {
        versionName = this.getCommitVersion(hash);
      } else {
        versionName = `h${index}-${hash}`;
      }
      
      try {
        const version = this.registry.addVersion(versionName, {
          description: message,
          auto: true,
          fromHistory: true,
          historyIndex: index,
          versionType: this.parseCommitMessage(message)
        });
        versions.push(version);
      } catch (error) {
        console.log(`⚠️ Skipped ${versionName}: ${error.message}`);
      }
    });
    
    return versions;
  }

  // Auto-add versions for major commits
  addMajorVersions(options = {}) {
    try {
      // Get commits with major changes
      const majorCommits = execSync(`git log --oneline --grep="feat\\|feature\\|breaking\\|major" -10`).toString().trim();
      
      if (!majorCommits) {
        console.log('No major commits found');
        return [];
      }
      
      const commits = majorCommits.split('\n');
      const versions = [];
      
      commits.forEach((commit, index) => {
        const [hash, ...messageParts] = commit.split(' ');
        const message = messageParts.join(' ');
        const versionName = options.useCommitNumber ? 
          this.getCommitVersion(hash) : 
          `major-${index}-${hash}`;
        
        try {
          const version = this.registry.addVersion(versionName, {
            description: message,
            auto: true,
            majorVersion: true,
            versionType: this.parseCommitMessage(message)
          });
          versions.push(version);
        } catch (error) {
          console.log(`⚠️ Skipped ${versionName}: ${error.message}`);
        }
      });
      
      return versions;
      
    } catch (error) {
      console.log('No major commits found or git error');
      return [];
    }
  }

  // Generate version name suggestions
  suggestVersionNames(commitHash = 'HEAD') {
    const commit = commitHash === 'HEAD' ? 
      execSync('git rev-parse HEAD').toString().trim() : commitHash;
    const message = execSync(`git log -1 --pretty=%B ${commit}`).toString().trim();
    
    return {
      commitBased: this.getCommitVersion(commit),
      semantic: this.getSemanticVersion(),
      combined: `${this.getSemanticVersion()}-${this.getCommitVersion(commit)}`,
      simple: commit.substring(0, 8),
      descriptive: this.generateDescriptiveName(message),
      timestamp: `v${new Date().toISOString().slice(0, 10).replace(/-/g, '')}`
    };
  }

  generateDescriptiveName(commitMessage) {
    // Extract meaningful words from commit message
    const words = commitMessage.toLowerCase()
      .replace(/[^\w\s]/g, ' ')
      .split(/\s+/)
      .filter(word => word.length > 2)
      .filter(word => !['and', 'the', 'for', 'with', 'from', 'add', 'fix', 'update'].includes(word))
      .slice(0, 3);
    
    return words.length > 0 ? words.join('-') : 'version';
  }
}

module.exports = AutoVersioning;
