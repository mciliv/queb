#!/usr/bin/env node
// Production build script - removes all debug code

const fs = require('fs-extra');
const path = require('path');
const { execSync } = require('child_process');

const BUILD_DIR = 'dist';
const DEBUG_PATTERNS = [
  /console\.(log|debug|warn|trace)\([^)]*\);?/g,
  /\/\*\s*DEBUG\s*START\s*\*\/[\s\S]*?\/\*\s*DEBUG\s*END\s*\*\//g,
  /\/\/\s*DEBUG:.*$/gm,
  /debug\.(log|warn|error|inspect|trace)\([^)]*\);?/g
];

const DEBUG_FILES = [
  'frontend/components/debug-events.js',
  'frontend/components/debug-utils.js',
  'backend/services/debug-utils.js',
  'test/debug-pipeline',
  'dev-debug'
];

async function cleanBuild() {
  console.log('🧹 Cleaning previous build...');
  await fs.remove(BUILD_DIR);
  await fs.ensureDir(BUILD_DIR);
}

async function copySource() {
  console.log('📂 Copying source files...');
  
  const excludePatterns = [
    'node_modules',
    '.git',
    'test',
    'tests',
    BUILD_DIR,
    ...DEBUG_FILES
  ];

  await fs.copy('.', BUILD_DIR, {
    filter: (src) => {
      const relativePath = path.relative('.', src);
      return !excludePatterns.some(pattern => relativePath.includes(pattern));
    }
  });
}

async function removeDebugCode() {
  console.log('🚫 Removing debug code...');
  
  const jsFiles = await getAllJSFiles(BUILD_DIR);
  
  for (const file of jsFiles) {
    let content = await fs.readFile(file, 'utf8');
    let originalContent = content;
    
    // Remove debug patterns
    DEBUG_PATTERNS.forEach(pattern => {
      content = content.replace(pattern, '');
    });
    
    // Remove debug imports
    content = content.replace(/import.*debug.*from.*;\s*/gi, '');
    content = content.replace(/require\(['"].*debug.*['"]\);\s*/gi, '');
    
    // Remove debug blocks
    content = content.replace(/if\s*\(\s*!?IS_PRODUCTION.*?\{[\s\S]*?\}/g, '');
    content = content.replace(/if\s*\(\s*DEBUG.*?\{[\s\S]*?\}/g, '');
    
    if (content !== originalContent) {
      await fs.writeFile(file, content);
      console.log(`   ✅ Cleaned: ${path.relative(BUILD_DIR, file)}`);
    }
  }
}

async function getAllJSFiles(dir) {
  const files = [];
  const entries = await fs.readdir(dir, { withFileTypes: true });
  
  for (const entry of entries) {
    const fullPath = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      files.push(...await getAllJSFiles(fullPath));
    } else if (entry.name.endsWith('.js')) {
      files.push(fullPath);
    }
  }
  
  return files;
}

async function createProductionEnv() {
  console.log('⚙️  Creating production environment...');
  
  const prodEnv = `NODE_ENV=production
DEBUG_MODE=false
CONSOLE_LOGS=false
`;
  
  await fs.writeFile(path.join(BUILD_DIR, '.env.production'), prodEnv);
}

async function installProdDependencies() {
  console.log('📦 Installing production dependencies...');
  
  const packagePath = path.join(BUILD_DIR, 'package.json');
  const packageJson = await fs.readJson(packagePath);
  
  // Remove dev dependencies
  delete packageJson.devDependencies;
  
  // Remove debug scripts
  const prodScripts = {};
  Object.entries(packageJson.scripts || {}).forEach(([key, value]) => {
    if (!key.includes('debug') && !key.includes('dev') && !key.includes('test')) {
      prodScripts[key] = value;
    }
  });
  packageJson.scripts = prodScripts;
  
  await fs.writeJson(packagePath, packageJson, { spaces: 2 });
  
  execSync('npm ci --production', { cwd: BUILD_DIR, stdio: 'inherit' });
}

async function main() {
  try {
    console.log('🚀 Starting production build...');
    
    await cleanBuild();
    await copySource();
    await removeDebugCode();
    await createProductionEnv();
    await installProdDependencies();
    
    console.log('✅ Production build complete!');
    console.log(`📁 Build output: ${BUILD_DIR}/`);
    console.log('🚫 All debug code removed');
    console.log('⚡ Ready for deployment');
    
  } catch (error) {
    console.error('❌ Build failed:', error.message);
    process.exit(1);
  }
}

if (require.main === module) {
  main();
}

module.exports = { main }; 