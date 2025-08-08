#!/usr/bin/env node

const esbuild = require('esbuild');
const path = require('path');
const fs = require('fs');

async function build() {
  const outdir = path.resolve(__dirname, '..', 'frontend', 'dist');
  const entry = path.resolve(__dirname, '..', 'frontend', 'core', 'index.jsx');
  fs.mkdirSync(outdir, { recursive: true });

  await esbuild.build({
    entryPoints: [entry],
    outfile: path.join(outdir, 'bundle.js'),
    bundle: true,
    sourcemap: true,
    minify: false,
    loader: { '.js': 'jsx', '.jsx': 'jsx' },
    define: { 'process.env.NODE_ENV': '"development"' },
  });

  console.log('✅ Frontend built to frontend/dist/bundle.js');
}

build().catch((err) => {
  console.error(err);
  process.exit(1);
});


