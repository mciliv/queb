#!/usr/bin/env node

// Simple SSL certificate generation script
const fs = require('fs');
const path = require('path');
const { execSync } = require('child_process');

const args = process.argv.slice(2);

// Default certificate directory
const certDir = args.find(arg => arg.startsWith('--dir=')) ?
  args.find(arg => arg.startsWith('--dir=')).replace('--dir=', '') :
  'config/certs';

console.log(`🔐 Generating SSL certificates in: ${certDir}`);

// Ensure directory exists
if (!fs.existsSync(certDir)) {
  fs.mkdirSync(certDir, { recursive: true });
  console.log(`📁 Created directory: ${certDir}`);
}

// Generate certificates using mkcert (assuming it's installed)
try {
  const certPath = path.join(certDir, 'localhost.pem');
  const keyPath = path.join(certDir, 'localhost-key.pem');

  console.log('📜 Generating localhost certificate...');
  execSync(`mkcert -cert-file "${certPath}" -key-file "${keyPath}" localhost 127.0.0.1 ::1`, {
    stdio: 'inherit'
  });

  console.log(`✅ Certificates generated:`);
  console.log(`   Certificate: ${certPath}`);
  console.log(`   Private Key: ${keyPath}`);

} catch (error) {
  console.error('❌ Failed to generate certificates. Make sure mkcert is installed:');
  console.error('   brew install mkcert');
  console.error('   mkcert -install');
  process.exit(1);
}
