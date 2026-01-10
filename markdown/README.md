# Queb - See the Molecules in Everything 🧪

AI-powered molecular analysis that reveals the chemical composition of anything you describe, photograph, or point your camera at.

## Features

- **Text Analysis**: Type any object → see its molecular components
- **Camera Analysis**: Point camera at objects → real-time molecular breakdown
- **Photo Analysis**: Upload images → discover chemical compositions
- **3D Visualization**: Interactive molecular structures (rotate, zoom, explore)

## Quick Start

### 1. Install Dependencies
```bash
npm install
```

### 2. Configure AI Providers
```bash
# Copy the example configuration
cp .env.example .env

# Edit .env and add your API keys:
# - OPENAI_API_KEY (from https://platform.openai.com/api-keys)
# - XAI_API_KEY (from https://console.x.ai/)
```

### 3. Validate Setup (Recommended)
```bash
# Validate AI configuration works across all environments
npm run validate:ai
```

### 4. Start Development Server
```bash
npm start              # http://localhost:8080
```

### 5. Switch AI Provider (Optional)
```bash
# Edit .env file
AI_PROVIDER=xai    # Switch to xAI/Grok
AI_PROVIDER=openai # Switch to OpenAI/GPT

# Restart the server
npm start
```

### 🔒 Locked-In Configuration

This setup is **environment-agnostic** and works consistently across:
- **Editors**: Cursor, Neovim, VSCode, any IDE¹
- **AI Providers**: OpenAI (GPT) or xAI (Grok)
- **Environments**: Local dev, production, testing

¹ *File ignoring is not synchronized between platforms (.cursorignore, .claudeignore, etc.)*

#### Environment Variables (`.env`)

```bash
# Choose your AI provider
AI_PROVIDER=openai    # or 'xai'

# OpenAI (uses 'latest' → automatically newest GPT model)
OPENAI_API_KEY=sk-...
OPENAI_MODEL=latest  # or specific: gpt-4, gpt-5

# xAI (uses 'latest' → automatically newest Grok model)
XAI_API_KEY=xai-...
XAI_MODEL=latest     # or specific: grok-3
```

#### Validation Command
```bash
npm run validate:ai   # Ensures everything works
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
