#!/usr/bin/env node
"use strict";

// Google AI Music CLI using Magenta.js (MusicVAE) - no fallbacks.
const fs = require('fs');
const path = require('path');

const mm = require('@magenta/music');  // Assume loads

async function main() {
  const args = process.argv.slice(2);
  const outArg = args.find(a => a.startsWith('--output='));
  const outputFile = path.resolve(process.cwd(), outArg ? outArg.split('=')[1] : 'output.mid');

  const modelUrl = 'https://storage.googleapis.com/magentadata/js/checkpoints/music_vae/mel_4bar_small';
  const MusicVAE = mm.MusicVAE;
  const model = new MusicVAE(modelUrl);
  await model.initialize();

  const [sequence] = await model.sample(1);
  const midiBytes = mm.sequenceProtoToMidi(sequence);
  fs.writeFileSync(outputFile, Buffer.from(midiBytes), 'binary');
  console.log(`Wrote MIDI to ${outputFile}`);
}

main().catch(err => {
  console.error('Error:', err.message);
  process.exit(1);
});


