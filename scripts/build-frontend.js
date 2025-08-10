#!/usr/bin/env node

const esbuild = require('esbuild');
const path = require('path');
const fs = require('fs');

async function build() {
  const outdir = path.resolve(__dirname, '..', 'frontend', 'dist');
  const entry = path.resolve(__dirname, '..', 'frontend', 'core', 'index.jsx');
  fs.mkdirSync(outdir, { recursive: true });

  const args = process.argv.slice(2);
  const watch = args.includes('--watch');

  if (watch) {
    const ctx = await esbuild.context({
      entryPoints: [entry],
      outfile: path.join(outdir, 'bundle.js'),
      bundle: true,
      sourcemap: true,
      minify: false,
      loader: { '.js': 'jsx', '.jsx': 'jsx' },
      define: { 'process.env.NODE_ENV': '"development"' },
    });
    await ctx.watch();
    console.log('👀 Frontend watch build started (esbuild)');
  } else {
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
}

build().catch((err) => {
  console.error(err);
  process.exit(1);
});


