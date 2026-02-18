#!/usr/bin/env node
"use strict";

// MIDI playback using Timidity live audio - no fallbacks.
const { spawnSync } = require('child_process');
const path = require('path');

function main() {
  const args = process.argv.slice(2);
  const midiPath = args[0] || 'output.mid';
  const resolved = path.resolve(process.cwd(), midiPath);

  // Direct Timidity live play
  const t = spawnSync('timidity', ['-iA', resolved], { stdio: 'inherit' });
  if (t.error || t.status !== 0) {
    console.error('Timidity failed - ensure installed and audio configured.');
    process.exit(1);
  }
  console.log('Played via Timidity.');
  process.exit(0);
}

main();


