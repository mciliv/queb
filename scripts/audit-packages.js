#!/usr/bin/env node

/**
 * Package Audit Script
 * Identifies unnecessary packages in package.json by scanning actual usage
 */

const fs = require('fs');
const path = require('path');
const { execSync } = require('child_process');

// Packages that are actually used in the codebase
const USED_PACKAGES = {
  // Core runtime
  'express': true,
  'cors': true,
  'dotenv': true,
  'openai': true,
  'zod': true,
  'pg': true,
  'livereload': true,
  'connect-livereload': true,
  
  // Testing
  'jest': true,
  'jest-environment-jsdom': true,
  'supertest': true,
  'jsdom': true,
  'puppeteer': true,
  'jest-fetch-mock': true,
  'chokidar': true,
  
  // Development
  'nodemon': true,
  'prettier': true,
  
  // Node.js built-ins (not in package.json)
  'fs': 'built-in',
  'path': 'built-in',
  'child_process': 'built-in',
  'http': 'built-in',
  'https': 'built-in',
  'net': 'built-in',
  'os': 'built-in',
  'util': 'built-in'
};

// Packages that can be safely removed
const UNNECESSARY_PACKAGES = [
  'multer',
  'formidable', 
  'busboy',
  'uuid',
  'ws',
  'prompts',
  'yargs',
  'yargs-parser',
  'body-parser'
];

function scanForPackageUsage(packageName) {
  const patterns = [
    `require('${packageName}')`,
    `require("${packageName}")`,
    `import.*from.*['"]${packageName}['"]`,
    `from ['"]${packageName}['"]`,
    `import ['"]${packageName}['"]`
  ];
  
  const directories = [
    'backend',
    'frontend', 
    'test',
    'scripts'
  ];
  
  let found = false;
  
  for (const dir of directories) {
    if (!fs.existsSync(dir)) continue;
    
    try {
      const result = execSync(`grep -r "${packageName}" ${dir} --include="*.js" --include="*.ts" --include="*.json" 2>/dev/null || true`, { encoding: 'utf8' });
      if (result.trim()) {
        // Check if it's actually a require/import, not just a comment or string
        const lines = result.split('\n').filter(line => line.trim());
        for (const line of lines) {
          if (patterns.some(pattern => line.includes(pattern.replace(/['"]/g, '').replace('require(', '').replace('import', '')))) {
            found = true;
            break;
          }
        }
      }
    } catch (error) {
      // Ignore grep errors
    }
  }
  
  return found;
}

function readPackageJson() {
  const packagePath = path.join(process.cwd(), 'package.json');
  return JSON.parse(fs.readFileSync(packagePath, 'utf8'));
}

function analyzePackageUsage() {
  const pkg = readPackageJson();
  const allDeps = { ...pkg.dependencies, ...pkg.devDependencies };
  const analysis = {
    used: [],
    unused: [],
    transitive: [],
    needsVerification: []
  };
  
  console.log('🔍 Scanning codebase for package usage...\n');
  
  for (const [name, version] of Object.entries(allDeps)) {
    if (USED_PACKAGES[name]) {
      analysis.used.push({ name, version, reason: 'Known essential' });
    } else if (UNNECESSARY_PACKAGES.includes(name)) {
      analysis.unused.push({ name, version, reason: 'Known unnecessary' });
    } else {
      // Check if it's actually used in the codebase
      const isUsed = scanForPackageUsage(name);
      if (isUsed) {
        analysis.used.push({ name, version, reason: 'Found in codebase' });
      } else {
        // Check if it's a Jest transitive dependency
        if (name.startsWith('jest-') || name.includes('istanbul') || name.includes('babel')) {
          analysis.transitive.push({ name, version, reason: 'Jest/Babel transitive' });
        } else {
          analysis.needsVerification.push({ name, version, reason: 'Needs verification' });
        }
      }
    }
  }
  
  return analysis;
}

function generateDetailedReport() {
  const analysis = analyzePackageUsage();
  
  console.log('📦 Detailed Package Audit Report\n');
  
  console.log('✅ ESSENTIAL PACKAGES:');
  analysis.used.forEach(pkg => {
    console.log(`  ✓ ${pkg.name}@${pkg.version} (${pkg.reason})`);
  });
  
  console.log('\n❌ UNNECESSARY PACKAGES:');
  analysis.unused.forEach(pkg => {
    console.log(`  ✗ ${pkg.name}@${pkg.version} (${pkg.reason})`);
  });
  
  console.log('\n🔄 TRANSITIVE DEPENDENCIES (Keep unless issues):');
  analysis.transitive.forEach(pkg => {
    console.log(`  🔄 ${pkg.name}@${pkg.version} (${pkg.reason})`);
  });
  
  if (analysis.needsVerification.length > 0) {
    console.log('\n⚠️  NEEDS VERIFICATION (Check before removing):');
    analysis.needsVerification.forEach(pkg => {
      console.log(`  ? ${pkg.name}@${pkg.version} (${pkg.reason})`);
    });
  }
  
  console.log('\n📊 STATISTICS:');
  console.log(`  Total packages: ${analysis.used.length + analysis.unused.length + analysis.transitive.length + analysis.needsVerification.length}`);
  console.log(`  Essential: ${analysis.used.length}`);
  console.log(`  Unnecessary: ${analysis.unused.length}`);
  console.log(`  Transitive: ${analysis.transitive.length}`);
  console.log(`  Needs verification: ${analysis.needsVerification.length}`);
  
  return analysis;
}

function suggestRemovals() {
  const analysis = generateDetailedReport();
  
  console.log('\n💡 REMOVAL SUGGESTIONS:');
  console.log('\n1. SAFE TO REMOVE (confirmed unused):');
  analysis.unused.forEach(pkg => {
    console.log(`   npm uninstall ${pkg.name}`);
  });
  
  if (analysis.needsVerification.length > 0) {
    console.log('\n2. VERIFY BEFORE REMOVING:');
    analysis.needsVerification.forEach(pkg => {
      console.log(`   # Check: ${pkg.name} - ${pkg.reason}`);
    });
  }
  
  console.log('\n3. KEEP (transitive dependencies):');
  analysis.transitive.forEach(pkg => {
    console.log(`   # Keep: ${pkg.name} - ${pkg.reason}`);
  });
}

function removeOnlyUnnecessary() {
  const analysis = analyzePackageUsage();
  const toRemove = analysis.unused.map(pkg => pkg.name);
  
  if (toRemove.length === 0) {
    console.log('✅ No unnecessary packages to remove');
    return;
  }
  
  console.log(`🗑️  Removing ${toRemove.length} confirmed unnecessary packages...`);
  console.log('Packages to remove:', toRemove.join(' '));
  
  try {
    execSync(`npm uninstall ${toRemove.join(' ')}`, { stdio: 'inherit' });
    console.log('✅ Successfully removed unnecessary packages');
  } catch (error) {
    console.error('❌ Failed to remove packages:', error.message);
  }
}

// CLI
const command = process.argv[2];

switch (command) {
  case 'report':
    generateDetailedReport();
    break;
  case 'suggest':
    suggestRemovals();
    break;
  case 'remove':
    removeOnlyUnnecessary();
    break;
  default:
    console.log('Usage: node scripts/audit-packages.js [report|suggest|remove]');
    console.log('  report  - Generate detailed package audit report');
    console.log('  suggest - Show removal suggestions');
    console.log('  remove  - Remove only confirmed unnecessary packages');
    break;
} 