#!/usr/bin/env node
"use strict";

// A minimal Google AI music CLI using Magenta.js (MusicVAE).
// Generates a short melody and writes it to output.mid.

const fs = require('fs');
const path = require('path');

async function main() {
  const args = process.argv.slice(2);
  // Basic argument parsing
  const modelUrlArg = args.find((a) => a.startsWith("--model="));
  const outputArg = args.find((a) => a.startsWith("--output="));
  const modelBaseUrl = modelUrlArg ? modelUrlArg.split("=")[1] :
    'https://storage.googleapis.com/magentadata/js/checkpoints/music_vae/mel_4bar_small';
  const outputFile = outputArg ? outputArg.split("=")[1] : 'output.mid';

  let mm;
  try {
    mm = require('@magenta/music');
  } catch (e) {
    console.error('Failed to load @magenta/music. Please install dependencies with:');
    console.error('  npm install @magenta/music @tensorflow/tfjs-node');
    process.exit(1);
  }

  try {
    // Initialize MusicVAE model
    const MusicVAE = mm.MusicVAE;
    if (!MusicVAE) {
      console.error('Unsupported Magenta MusicVAE entry. Ensure correct version of @magenta/music is installed.');
      process.exit(1);
    }

    const model = new MusicVAE(modelBaseUrl);
    await model.initialize();

    // Generate a single melody
    const samples = await model.sample(1);
    const sequenceProto = samples[0];

    // Convert to MIDI bytes
    const midiBytes = mm.sequenceProtoToMidi(sequenceProto);
    const outputPath = path.resolve(process.cwd(), outputFile);
    fs.writeFileSync(outputPath, Buffer.from(midiBytes), 'binary');
    console.log(`Wrote MIDI to ${outputPath}`);
    process.exit(0);
  } catch (err) {
    console.error('Error generating music:', err && err.message ? err.message : err);
    process.exit(1);
  }
}

main();







