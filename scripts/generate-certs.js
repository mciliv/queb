#!/usr/bin/env node

// Wrapper script to generate SSL certificates using the generic toolkit
const { spawn } = require('child_process');
const path = require('path');

const args = process.argv.slice(2);

// Use config/certs as the default directory
const certDir = args.find(arg => arg.startsWith('--dir=')) ?
  args.find(arg => arg.startsWith('--dir=')) :
  '--dir=config/certs';

const filteredArgs = args.filter(arg => !arg.startsWith('--dir='));
const finalArgs = [certDir, ...filteredArgs];

// Find util directory (assume it's at same level as mol)
const utilPath = path.dirname(path.dirname(__dirname));
const toolkitPath = path.join(utilPath, 'util', 'dev-toolkit', 'certs', 'generate-certs.js');

// Execute the generic certificate generator
const child = spawn('node', [toolkitPath, ...finalArgs], {
  stdio: 'inherit',
  cwd: path.dirname(__dirname) // Set cwd to mol project root
});

child.on('exit', (code) => {
  process.exit(code);
});
