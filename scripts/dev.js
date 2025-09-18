#!/usr/bin/env node

// Simple development server script
const { spawn } = require('child_process');
const path = require('path');

console.log('🚀 Starting Molecular Analysis Development Server...');

// Start the server directly
const serverPath = path.join(__dirname, '..', 'index.js');
const child = spawn('node', [serverPath], {
  stdio: 'inherit',
  cwd: path.join(__dirname, '..')
});

child.on('exit', (code) => {
  if (code === 0) {
    console.log('✅ Development server exited cleanly');
  } else {
    console.log(`❌ Development server exited with code ${code}`);
  }
  process.exit(code);
});

child.on('error', (error) => {
  console.error('❌ Failed to start development server:', error.message);
  process.exit(1);
});

// Handle parent process signals
process.on('SIGINT', () => {
  console.log('\n🛑 Received SIGINT, shutting down server...');
  child.kill('SIGINT');
});

process.on('SIGTERM', () => {
  console.log('\n🛑 Received SIGTERM, shutting down server...');
  child.kill('SIGTERM');
});






