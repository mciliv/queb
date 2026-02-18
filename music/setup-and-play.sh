#!/bin/bash

# Setup and play script for music-playground
# Installs dependencies if needed and runs generation + playback

set -e

echo "=== Music Playground Setup & Play ==="

# Ensure we're in the right dir
cd "$(dirname "$0")"

# Install timidity if not present (use absolute brew path)
if ! command -v timidity &> /dev/null; then
    echo "Installing timidity via Homebrew..."
    /opt/homebrew/bin/brew install timidity || {
        echo "Failed to install timidity. Continuing without CLI playback."
    }
fi

# Install npm deps if not done
if [ ! -d "node_modules" ]; then
    echo "Installing Node dependencies..."
    npm install
fi

# Generate MIDI
echo "Generating MIDI..."
npm run generate

# Play MIDI
echo "Playing MIDI..."
npm run play

echo "Done! Check output.mid if playback didn't work."





