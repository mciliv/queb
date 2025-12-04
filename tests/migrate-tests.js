#!/usr/bin/env node

/**
 * Migration script to move existing tests to the new structure
 * This script will:
 * 1. Move test files to appropriate __tests__ subdirectories
 * 2. Update import paths
 * 3. Rename files to follow new naming conventions
 * 4. Create a migration report
 */

const fs = require('fs');
const path = require('path');

// Test file mappings
const FILE_MAPPINGS = {
  // Unit tests
  'suites/unit/*.test.js': 'unit',
  'unit/*.test.js': 'unit',
  
  // Integration tests
  'suites/integration/*.test.js': 'integration',
  'integration/*.test.js': 'integration',
  
  // E2E tests
  'suites/e2e/*.test.js': 'e2e',
  
  // Visual tests
  'suites/visual/*.test.js': 'visual',
  'visual/*.test.js': 'visual'
};

// Files to skip
const SKIP_FILES = [
  'jest.config.js',
  'migrate-tests.js',
  'TEST_REFACTORING_PLAN.md',
  'README.md'
];

class TestMigrator {
  constructor() {
    this.report = {
      migrated: [],
      skipped: [],
      errors: [],
      warnings: []
    };
  }

  async migrate() {
    console.log('Starting test migration...\n');
    
    // Create backup
    await this.createBackup();
    
    // Find all test files
    const testFiles = this.findTestFiles();
    console.log(`Found ${testFiles.length} test files to process\n`);
    
    // Process each file
    for (const file of testFiles) {
      await this.processFile(file);
    }
    
    // Generate report
    this.generateReport();
  }

  createBackup() {
    const backupDir = path.join(__dirname, '__backup__', new Date().toISOString().replace(/:/g, '-'));
    
    console.log(`Creating backup in ${backupDir}...`);
    
    // Copy current test directory structure
    this.copyDirectory(__dirname, backupDir, [
      'node_modules',
      '__backup__',
      'coverage',
      'reports',
      '__tests__',
      '__config__',
      '__utils__',
      '__fixtures__'
    ]);
    
    console.log('Backup created successfully\n');
    
    this.report.backupLocation = backupDir;
  }

  copyDirectory(src, dest, excludes = []) {
    if (!fs.existsSync(dest)) {
      fs.mkdirSync(dest, { recursive: true });
    }
    
    const entries = fs.readdirSync(src, { withFileTypes: true });
    
    for (const entry of entries) {
      const srcPath = path.join(src, entry.name);
      const destPath = path.join(dest, entry.name);
      
      if (excludes.includes(entry.name)) continue;
      
      if (entry.isDirectory()) {
        this.copyDirectory(srcPath, destPath, excludes);
      } else {
        fs.copyFileSync(srcPath, destPath);
      }
    }
  }

  findTestFiles() {
    const files = [];
    
    const walk = (dir) => {
      const entries = fs.readdirSync(dir, { withFileTypes: true });
      
      for (const entry of entries) {
        const fullPath = path.join(dir, entry.name);
        const relativePath = path.relative(__dirname, fullPath);
        
        if (entry.isDirectory()) {
          // Skip new structure directories
          if (!['__tests__', '__config__', '__utils__', '__fixtures__', '__backup__', 'node_modules'].includes(entry.name)) {
            walk(fullPath);
          }
        } else if (entry.isFile() && entry.name.endsWith('.test.js')) {
          if (!SKIP_FILES.includes(entry.name)) {
            files.push(relativePath);
          }
        }
      }
    };
    
    walk(__dirname);
    return files;
  }

  async processFile(filePath) {
    console.log(`Processing: ${filePath}`);
    
    try {
      // Determine target directory
      const targetDir = this.determineTargetDir(filePath);
      if (!targetDir) {
        this.report.skipped.push({
          file: filePath,
          reason: 'Could not determine test type'
        });
        console.log(`  ⚠️  Skipped - could not determine test type\n`);
        return;
      }
      
      // Read file content
      const content = fs.readFileSync(path.join(__dirname, filePath), 'utf8');
      
      // Update imports
      const updatedContent = this.updateImports(content, filePath);
      
      // Determine new filename
      const newFilename = this.getNewFilename(path.basename(filePath));
      const newPath = path.join(__dirname, '__tests__', targetDir, newFilename);
      
      // Ensure directory exists
      const newDir = path.dirname(newPath);
      if (!fs.existsSync(newDir)) {
        fs.mkdirSync(newDir, { recursive: true });
      }
      
      // Write updated file
      fs.writeFileSync(newPath, updatedContent);
      
      this.report.migrated.push({
        from: filePath,
        to: path.relative(__dirname, newPath),
        type: targetDir
      });
      
      console.log(`  ✅ Migrated to __tests__/${targetDir}/${newFilename}\n`);
      
    } catch (error) {
      this.report.errors.push({
        file: filePath,
        error: error.message
      });
      console.log(`  ❌ Error: ${error.message}\n`);
    }
  }

  determineTargetDir(filePath) {
    // Check explicit mappings
    for (const [pattern, type] of Object.entries(FILE_MAPPINGS)) {
      const regex = new RegExp(pattern.replace('*', '.*'));
      if (regex.test(filePath)) {
        return type;
      }
    }
    
    // Try to determine from content
    const content = fs.readFileSync(path.join(__dirname, filePath), 'utf8');
    
    if (content.includes('puppeteer') || content.includes('browser.newPage')) {
      return 'e2e';
    }
    if (content.includes('screenshot') || content.includes('visual')) {
      return 'visual';
    }
    if (content.includes('supertest') || content.includes('integration')) {
      return 'integration';
    }
    
    // Default to unit
    return 'unit';
  }

  getNewFilename(filename) {
    // Convert *.test.js to *.spec.js
    return filename.replace('.test.js', '.spec.js');
  }

  updateImports(content, filePath) {
    let updated = content;
    
    // Update relative imports to use path aliases
    const replacements = [
      // Source imports
      [/require\(['"]\.\.\/\.\.\/\.\.\/src\//g, "require('@/"],
      [/require\(['"]\.\.\/\.\.\/src\//g, "require('@/"],
      [/from ['"]\.\.\/\.\.\/\.\.\/src\//g, "from '@/"],
      [/from ['"]\.\.\/\.\.\/src\//g, "from '@/"],
      
      // Test utility imports
      [/require\(['"]\.\.\/\.\.\/support\/fixtures\//g, "require('@fixtures/"],
      [/require\(['"]\.\.\/support\/fixtures\//g, "require('@fixtures/"],
      [/from ['"]\.\.\/\.\.\/support\/fixtures\//g, "from '@fixtures/"],
      [/from ['"]\.\.\/support\/fixtures\//g, "from '@fixtures/"],
      
      // Update to new test utilities
      [/require\(['"]\.\.\/\.\.\/support\/fixtures\/utils['"]\)/g, "require('@test/helpers')"],
      [/from ['"]\.\.\/\.\.\/support\/fixtures\/utils['"]/g, "from '@test/helpers'"],
    ];
    
    for (const [pattern, replacement] of replacements) {
      updated = updated.replace(pattern, replacement);
    }
    
    // Add custom matcher import if using expect extensions
    if (updated.includes('expect.extend') || updated.includes('toBeValidSmiles')) {
      if (!updated.includes('@test/matchers')) {
        const lines = updated.split('\n');
        const lastRequireIndex = lines.findLastIndex(line => 
          line.includes('require(') || line.includes('import ')
        );
        
        if (lastRequireIndex !== -1) {
          lines.splice(lastRequireIndex + 1, 0, 
            "const { extendExpect } = require('@test/matchers');\n\n// Extend Jest with custom matchers\nextendExpect();"
          );
          updated = lines.join('\n');
        }
      }
    }
    
    return updated;
  }

  generateReport() {
    console.log('\n' + '='.repeat(60));
    console.log('MIGRATION REPORT');
    console.log('='.repeat(60) + '\n');
    
    console.log(`Backup Location: ${this.report.backupLocation}\n`);
    
    console.log(`Successfully Migrated: ${this.report.migrated.length}`);
    this.report.migrated.forEach(m => {
      console.log(`  ${m.from} → ${m.to} (${m.type})`);
    });
    
    if (this.report.skipped.length > 0) {
      console.log(`\nSkipped: ${this.report.skipped.length}`);
      this.report.skipped.forEach(s => {
        console.log(`  ${s.file} - ${s.reason}`);
      });
    }
    
    if (this.report.errors.length > 0) {
      console.log(`\nErrors: ${this.report.errors.length}`);
      this.report.errors.forEach(e => {
        console.log(`  ${e.file} - ${e.error}`);
      });
    }
    
    if (this.report.warnings.length > 0) {
      console.log(`\nWarnings: ${this.report.warnings.length}`);
      this.report.warnings.forEach(w => {
        console.log(`  ${w.file} - ${w.warning}`);
      });
    }
    
    // Save report
    const reportPath = path.join(__dirname, 'migration-report.json');
    fs.writeFileSync(reportPath, JSON.stringify(this.report, null, 2));
    console.log(`\nDetailed report saved to: ${reportPath}`);
    
    console.log('\n' + '='.repeat(60));
    console.log('NEXT STEPS:');
    console.log('='.repeat(60));
    console.log('1. Review migrated files for any manual adjustments needed');
    console.log('2. Run tests to ensure they still pass: npm test');
    console.log('3. Update any CI/CD configurations');
    console.log('4. Delete old test files once migration is verified');
    console.log('5. Remove the backup directory when no longer needed\n');
  }
}

// Run migration if called directly
if (require.main === module) {
  const migrator = new TestMigrator();
  
  console.log('Test Migration Tool\n');
  console.log('This will migrate your existing tests to the new structure.');
  console.log('A backup will be created before any changes are made.\n');
  
  const readline = require('readline').createInterface({
    input: process.stdin,
    output: process.stdout
  });
  
  readline.question('Do you want to proceed? (y/N) ', (answer) => {
    readline.close();
    
    if (answer.toLowerCase() === 'y') {
      migrator.migrate().catch(error => {
        console.error('Migration failed:', error);
        process.exit(1);
      });
    } else {
      console.log('Migration cancelled');
    }
  });
}

module.exports = TestMigrator;
