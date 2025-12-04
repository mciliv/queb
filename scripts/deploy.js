#!/usr/bin/env node

const { spawn } = require('child_process');
const path = require('path');

function findProjectRoot() {
  const scriptDir = __dirname;
  return path.dirname(scriptDir);
}

async function deploy(projectRoot, env = 'dev') {
  console.log(`🚀 Deploying to ${env === 'dev' ? 'dev.queb.space' : 'queb.space'}...`);
  
  const deployScript = path.join(projectRoot, 'deploy', 'deploy.sh');
  
  return new Promise((resolve, reject) => {
    const child = spawn('bash', [deployScript, env], {
      stdio: 'inherit',
      cwd: projectRoot,
      env: process.env
    });
    
    child.on('exit', (code) => {
      if (code === 0) {
        resolve();
      } else {
        reject(new Error(`Deploy script exited with code ${code}`));
      }
    });
    
    child.on('error', reject);
  });
}

async function main() {
  const projectRoot = findProjectRoot();
  const env = process.argv[2] || 'dev';
  
  if (env !== 'dev' && env !== 'prod') {
    console.error('Usage: node scripts/deploy.js [dev|prod]');
    console.error('  dev  - Deploy to dev.queb.space');
    console.error('  prod - Deploy to queb.space');
    process.exit(1);
  }
  
  await deploy(projectRoot, env);
}

main().catch((error) => {
  console.error('❌ Deployment failed:', error.message);
  process.exit(1);
});
