# Music Playground

A standalone repository to generate and play MIDI music produced by Google AI Music tools.

What’s included
- A minimal Node.js project with a MusicVAE playback workflow (`music-playground/play.js`).
- A CLI-friendly package.json script to play the generated MIDI via Timidity when available.
- A simple README with usage instructions.

Prerequisites
- Node.js >= 14
- Optional: Timidity (for MIDI playback) installed on the host system and available on PATH
- If you want to generate music with Magenta, you’ll also need @magenta/music and @tensorflow/tfjs-node installed in this repo when you extend functionality

Getting started
1) Install dependencies
```
cd music-playground
npm install
```

2) Generate a MIDI file (optional: you can use the existing script to simply generate an output.mid)
- The repository currently ships with a basic playback script; you can modify or extend it to call into Magenta as needed in your environment.

3) Play the MIDI file (requires Timidity or your preferred MIDI player)
```
npm run play
```

Notes
- The generated MIDI is written to `output.mid` by default. You can customize the file name using the CLI (see `play.js`).
- The `.gitignore` excludes `output.mid` by default to avoid committing large binary files.







