#!/usr/bin/env node

const esbuild = require('esbuild');
const path = require('path');
const fs = require('fs');

async function build() {
  // Output to frontend/dist so server serves /dist/* correctly
  const outdir = path.resolve(__dirname, 'dist');
  const entry = path.resolve(__dirname, 'core', 'index.jsx');
  fs.mkdirSync(outdir, { recursive: true });

  const args = process.argv.slice(2);
  const watch = args.includes('--watch');
  const quiet = args.includes('--quiet') || args.includes('-q');
  const dev = args.includes('--dev') || process.env.NODE_ENV === 'development';

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
      logLevel: quiet ? 'silent' : 'info',
      external: [
        'config/*', 
        '../../config/*',
        'config/project.js',
        'config/env.js',
        'backend/*',
        'src/server/*'
      ],
      banner: {
        js: `
          // Development mode optimizations
          if (typeof window !== 'undefined') {
            console.log('[BUILD] Frontend bundle loaded at:', new Date().toISOString());
          }
        `
      }
    });
    // Ensure initial build so /dist/bundle.js exists
    await ctx.rebuild();

    // Watch for changes - Chrome auto-reload will handle refreshing
    await ctx.watch(() => {
      if (!quiet) console.log('🔄 Frontend rebuilt');
    });

    if (!quiet) console.log('👀 Frontend watch build started (esbuild)');
  } else {
    await esbuild.build({
      entryPoints: [entry],
      outfile: path.join(outdir, 'bundle.js'),
      bundle: true,
      sourcemap: dev,
      minify: !dev,
      treeShaking: !dev,
      loader: { '.js': 'jsx', '.jsx': 'jsx' },
      define: {
        'process.env.NODE_ENV': dev ? '"development"' : '"production"',
        'process.env.REACT_APP_RUN_VISUAL_TESTS': '"false"',
        'global': 'globalThis'
      },
      platform: 'browser',
      drop: dev ? [] : ['console', 'debugger'],
      logLevel: quiet ? 'silent' : 'info',
      external: [
        'config/*', 
        '../../config/*',
        'config/project.js',
        'config/env.js',
        'backend/*',
        'src/server/*'
      ]
    });
  
    if (!quiet) console.log('✅ Frontend built to dist/bundle.js');
  }
}

build().catch((err) => {
  console.error(err);
  process.exit(1);
});
