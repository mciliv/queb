const VersionRegistry = require('../support/version-registry');
const VersionManager = require('./version-manager');
const fs = require('fs');
const path = require('path');

class VersionComparator {
  constructor() {
    this.registry = new VersionRegistry();
    this.manager = new VersionManager();
    this.resultsDir = path.join(__dirname, '../results');
    this.ensureResultsDir();
  }

  ensureResultsDir() {
    if (!fs.existsSync(this.resultsDir)) {
      fs.mkdirSync(this.resultsDir, { recursive: true });
    }
  }

  // Run comparison between two versions
  async compareVersions(version1, version2, testCases = ['water', 'ethanol', 'caffeine']) {
    console.log(`🔬 Comparing ${version1} vs ${version2}`);
    
    const results = {
      version1: { name: version1, results: {} },
      version2: { name: version2, results: {} },
      comparison: {},
      timestamp: new Date().toISOString()
    };

    try {
      // Start both versions
      await this.manager.startVersion(version1);
      await this.manager.startVersion(version2);
      
      // Connect browsers
      const browser1 = await this.manager.connectBrowser(version1);
      const browser2 = await this.manager.connectBrowser(version2);
      
      // Run tests on both versions in parallel
      const testPromises = testCases.map(async (testCase) => {
        const [result1, result2] = await Promise.all([
          this.runSingleTest(browser1, version1, testCase),
          this.runSingleTest(browser2, version2, testCase)
        ]);
        
        results.version1.results[testCase] = result1;
        results.version2.results[testCase] = result2;
        results.comparison[testCase] = this.diffResults(result1, result2);
      });

      await Promise.all(testPromises);
      
      // Generate comparison report
      const reportPath = await this.generateComparisonReport(results);
      console.log(`📊 Comparison report: ${reportPath}`);
      
      return results;
      
    } catch (error) {
      console.error(`❌ Comparison failed: ${error.message}`);
      throw error;
    }
  }

  async runSingleTest(browser, version, testCase) {
    const config = this.registry.getVersion(version);
    const page = await browser.newPage();
    
    try {
      await page.goto(config.ports.http ? `http://localhost:${config.ports.http}` : 'http://localhost:3000', 
        { waitUntil: 'networkidle0' });
      
      // Wait for app to load
      await page.waitForSelector('input[type="text"]', { timeout: 10000 });
      
      // Input test case
      await page.type('input[type="text"]', testCase, { delay: 50 });
      await page.keyboard.press('Enter');
      
      // Wait for results
      await page.waitForSelector('.molecule-container', { timeout: 15000 });
      
      // Extract results
      const result = await page.evaluate(() => {
        const containers = document.querySelectorAll('.molecule-container');
        const molecules = [];
        
        containers.forEach(container => {
          const name = container.querySelector('.molecule-name')?.textContent || 'Unknown';
          const formula = container.querySelector('.molecule-formula')?.textContent || 'Unknown';
          molecules.push({ name, formula });
        });
        
        return {
          moleculeCount: containers.length,
          molecules,
          timestamp: new Date().toISOString()
        };
      });
      
      // Take screenshot
      const screenshotPath = path.join(this.resultsDir, `${version}_${testCase}_${Date.now()}.png`);
      await page.screenshot({ path: screenshotPath, fullPage: true });
      result.screenshot = screenshotPath;
      
      return result;
      
    } finally {
      await page.close();
    }
  }

  diffResults(result1, result2) {
    const diff = {
      moleculeCountDiff: result2.moleculeCount - result1.moleculeCount,
      moleculesAdded: [],
      moleculesRemoved: [],
      moleculesChanged: []
    };

    // Simple molecule comparison by name
    const molecules1 = result1.molecules.map(m => m.name);
    const molecules2 = result2.molecules.map(m => m.name);
    
    diff.moleculesAdded = molecules2.filter(m => !molecules1.includes(m));
    diff.moleculesRemoved = molecules1.filter(m => !molecules2.includes(m));
    
    // Check for formula changes
    result1.molecules.forEach(mol1 => {
      const mol2 = result2.molecules.find(m => m.name === mol1.name);
      if (mol2 && mol1.formula !== mol2.formula) {
        diff.moleculesChanged.push({
          name: mol1.name,
          oldFormula: mol1.formula,
          newFormula: mol2.formula
        });
      }
    });

    return diff;
  }

  async generateComparisonReport(results) {
    const reportPath = path.join(this.resultsDir, `comparison_${results.version1.name}_vs_${results.version2.name}_${Date.now()}.json`);
    
    // Enhanced report with analysis
    const report = {
      ...results,
      summary: this.generateSummary(results),
      recommendations: this.generateRecommendations(results)
    };
    
    fs.writeFileSync(reportPath, JSON.stringify(report, null, 2));
    
    // Also generate human-readable report
    const readableReport = this.generateReadableReport(report);
    const readablePath = reportPath.replace('.json', '.md');
    fs.writeFileSync(readablePath, readableReport);
    
    return readablePath;
  }

  generateSummary(results) {
    const testCases = Object.keys(results.comparison);
    let totalDiffs = 0;
    let improvements = 0;
    let regressions = 0;

    testCases.forEach(testCase => {
      const diff = results.comparison[testCase];
      if (diff.moleculeCountDiff > 0) improvements++;
      if (diff.moleculeCountDiff < 0) regressions++;
      totalDiffs += Math.abs(diff.moleculeCountDiff);
    });

    return {
      totalTestCases: testCases.length,
      totalDifferences: totalDiffs,
      improvements,
      regressions,
      overallTrend: improvements > regressions ? 'improvement' : 
                   regressions > improvements ? 'regression' : 'neutral'
    };
  }

  generateRecommendations(results) {
    const recommendations = [];
    const summary = results.summary;

    if (summary.regressions > 0) {
      recommendations.push('⚠️ Investigate regressions in molecular detection');
    }
    if (summary.improvements > 0) {
      recommendations.push('✅ Improvements detected - consider promoting version');
    }
    if (summary.totalDifferences === 0) {
      recommendations.push('🔄 No differences found - versions are equivalent');
    }

    return recommendations;
  }

  generateReadableReport(report) {
    const { version1, version2, comparison, summary } = report;
    
    let markdown = `# Version Comparison Report\n\n`;
    markdown += `**${version1.name}** vs **${version2.name}**\n`;
    markdown += `Generated: ${report.timestamp}\n\n`;
    
    markdown += `## Summary\n`;
    markdown += `- Test cases: ${summary.totalTestCases}\n`;
    markdown += `- Total differences: ${summary.totalDifferences}\n`;
    markdown += `- Improvements: ${summary.improvements}\n`;
    markdown += `- Regressions: ${summary.regressions}\n`;
    markdown += `- Overall trend: ${summary.overallTrend}\n\n`;

    markdown += `## Test Results\n\n`;
    
    Object.entries(comparison).forEach(([testCase, diff]) => {
      markdown += `### ${testCase}\n`;
      markdown += `- Molecule count change: ${diff.moleculeCountDiff > 0 ? '+' : ''}${diff.moleculeCountDiff}\n`;
      
      if (diff.moleculesAdded.length > 0) {
        markdown += `- Added: ${diff.moleculesAdded.join(', ')}\n`;
      }
      if (diff.moleculesRemoved.length > 0) {
        markdown += `- Removed: ${diff.moleculesRemoved.join(', ')}\n`;
      }
      if (diff.moleculesChanged.length > 0) {
        markdown += `- Changed formulas: ${diff.moleculesChanged.map(m => `${m.name} (${m.oldFormula} → ${m.newFormula})`).join(', ')}\n`;
      }
      markdown += `\n`;
    });

    markdown += `## Recommendations\n\n`;
    report.recommendations.forEach(rec => {
      markdown += `- ${rec}\n`;
    });

    return markdown;
  }

  // List all versions with their status
  listVersions() {
    const registry = this.loadRegistry();
    return Object.values(registry.versions).sort((a, b) => 
      new Date(b.timestamp) - new Date(a.timestamp)
    );
  }

  // Set current version for testing
  setCurrentVersion(versionName) {
    const registry = this.loadRegistry();
    if (!registry.versions[versionName]) {
      throw new Error(`Version ${versionName} not found`);
    }
    registry.currentVersion = versionName;
    this.saveRegistry(registry);
  }

  getCurrentVersion() {
    const registry = this.loadRegistry();
    return registry.currentVersion;
  }
}

module.exports = VersionComparator;
