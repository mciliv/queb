#!/usr/bin/env node

const esbuild = require('esbuild');
const path = require('path');
const fs = require('fs');

async function build() {
  const outdir = path.resolve(__dirname, '..', 'dist');
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
      define: { 
        'process.env.NODE_ENV': '"development"',
        'process.env.REACT_APP_RUN_VISUAL_TESTS': '"false"',
        'global': 'globalThis'
      },
      platform: 'browser',
    });
    await ctx.watch();
    console.log('👀 Fron‘tend watch build started (esbuild)');
  } else {
    await esbuild.build({
      entryPoints: [entry],
      outfile: path.join(outdir, 'bundle.js'),
      bundle: true,
      sourcemap: false,
      minify: true,
      treeShaking: true,
      loader: { '.js': 'jsx', '.jsx': 'jsx' },
      define: { 
        'process.env.NODE_ENV': '"production"',
        'process.env.REACT_APP_RUN_VISUAL_TESTS': '"false"',
        'global': 'globalThis'
      },
      platform: 'browser',
      drop: ['console', 'debugger'],
    });
  
    console.log('✅ Frontend built to dist/bundle.js');
  }
}

build().catch((err) => {
  console.error(err);
  process.exit(1);
});
