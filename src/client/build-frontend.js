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
      external: ['config/*', '../../config/*'],
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

    // Watch for changes - manual refresh required
    await ctx.watch(() => {
      if (!quiet) console.log('🔄 Frontend rebuilt - refresh browser to see changes');
    });

    if (!quiet) console.log('👀 Frontend watch build started (esbuild)');
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
      logLevel: quiet ? 'silent' : 'info',
      external: ['config/*', '../../config/*']
    });
  
    if (!quiet) console.log('✅ Frontend built to dist/bundle.js');
  }
}

build().catch((err) => {
  console.error(err);
  process.exit(1);
});
