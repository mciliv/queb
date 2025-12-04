# Queb - See the Molecules in Everything 🧪

AI-powered molecular analysis that reveals the chemical composition of anything you describe, photograph, or point your camera at.

## Features

- **Text Analysis**: Type any object → see its molecular components
- **Camera Analysis**: Point camera at objects → real-time molecular breakdown
- **Photo Analysis**: Upload images → discover chemical compositions
- **3D Visualization**: Interactive molecular structures (rotate, zoom, explore)

## Quick Start

```
npm install
cp .env.example .env   # add OPENAI_API_KEY
npm start              # http://localhost:8080
```

See `npm run` for all available scripts.

## Architecture

```
User Input → AI Analysis → Chemical Recognition → 
Structure Resolution → 3D Generation → Interactive Display
```

| Layer | Technology |
|-------|------------|
| Frontend | React 19, 3Dmol.js |
| Backend | Node.js, Express |
| AI | OpenAI GPT-4o |
| Chemistry | PubChem API, RDKit |
| Database | PostgreSQL (optional) |

## Project Structure

```
src/
├── client/      # React frontend
│   ├── components/  # UI components
│   ├── hooks/       # Custom React hooks
│   └── utils/       # Frontend utilities
├── server/      # Express backend
│   ├── api/         # API endpoints
│   └── services/    # Business logic
└── core/        # Shared modules
```

## How It Works

1. User provides text/image/camera input
2. GPT-4o identifies objects and their chemical components
3. Names resolved to SMILES via PubChem
4. SMILES converted to 3D SDF files using RDKit
5. 3Dmol.js renders interactive structures

## Configuration

Required in `.env`:
- `OPENAI_API_KEY` - OpenAI API key

Advanced settings in `config/config.toml`.

## Troubleshooting

| Issue | Fix |
|-------|-----|
| "No molecules found" | Be more specific (e.g., "caffeine" not "stimulant") |
| 3D viewer not loading | Enable WebGL, refresh, check console |
| Camera not working | Grant permissions, use HTTPS |
