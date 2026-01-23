# Chemical Analyzer

Display chemical structures contained in objects you specify or identify.

> **Note**: This repository contains all relevant projects in a monorepo structure.

## What It Does

Upload an image, take a photo, enter a link, or describe an object → Get AI-identified molecules → See interactive 3D molecular structures

**Example**: Photo of coffee → Identifies caffeine, water, etc. → Displays 3D molecular models

## Quick Start

1. **Install dependencies:**
   ```bash
   npm install
   ```

2. **Set up your API key:**
   Copy `.env.example` to `.env` and add your OpenAI API key:
   ```bash
   cp .env.example .env
   # Edit .env and add: OPENAI_API_KEY=sk-...
   ```

3. **Run the app:**
   ```bash
   npm start
   ```

4. **Open your browser:**
   Visit `http://localhost:8080`

## Input Methods

The app supports four ways to analyze objects:

- **Text**: Type object name (e.g., "coffee", "aspirin", "water")
- **Camera**: Point camera at objects for real-time analysis
- **Photo**: Upload images to analyze
- **Link**: Enter image URL to analyze

All methods identify chemical compounds and display them as interactive 3D molecular structures.

## How It Works

1. **Input**: You provide text, image, photo, or URL
2. **AI Analysis**: AI identifies the object and its chemical components
3. **Chemical Lookup**: Names are resolved to molecular structures via PubChem
4. **3D Generation**: Structures are converted to 3D coordinates
5. **Visualization**: 3Dmol.js renders interactive molecular structures (rotate, zoom, explore)

## API Endpoints

### POST /api/structuralize
Analyzes text for chemical compounds.

**Request:**
```json
{
  "text": "coffee",
  "lookupMode": "GPT-5"
}
```

**Response:**
```json
{
  "object": "coffee",
  "chemicals": [
    {
      "name": "Caffeine",
      "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    }
  ]
}
```

### POST /api/structuralize-image
Analyzes image for chemical compounds.

**Request:**
```json
{
  "imageBase64": "base64-encoded-image",
  "x": 150,
  "y": 200
}
```

### POST /api/generate-sdfs
Converts SMILES notation to SDF files for 3D display.

**Request:**
```json
{
  "smiles": ["CN1C=NC2=C1C(=O)N(C(=O)N2C)C"]
}
```

## Project Structure

This monorepo contains all relevant projects:

```
/
├── server.js              # Main server entry point
├── package.json           # Dependencies
├── README.md              # This file
├── .env.example          # Environment variable template
├── public/               # Static files
│   ├── manifest.json
│   └── sw.js
├── mcp/                  # Model Context Protocol server
│   └── chemical-discovery-server.js  # MCP server for chemical discovery
├── src/
│   ├── client/           # Frontend React app
│   │   ├── components/   # UI components
│   │   ├── hooks/        # React hooks
│   │   └── assets/       # CSS, icons
│   └── server/           # Backend services
│       ├── services/     # Core business logic
│       └── routes/       # API routes
└── resources/            # Data files
    └── chemical-databases.json
```

## Technology Stack

- **Frontend**: React, 3Dmol.js
- **Backend**: Node.js, Express
- **AI**: OpenAI API (GPT-4 Vision)
- **Chemistry**: PubChem API, SMILES notation, SDF files

## Development

This is a simplified, focused application for analyzing chemical contents and visualizing molecular structures.

Key features:
- Single Express server file
- React frontend with 3D visualization
- Direct AI API integration
- No database dependencies
- Simple, clear architecture

## Environment Variables

Create a `.env` file with:

```bash
OPENAI_API_KEY=sk-...your-key-here...
AI_PROVIDER=openai
PORT=8080
```

For xAI support:
```bash
AI_PROVIDER=xai
XAI_API_KEY=xai-...your-key-here...
```

## License

MIT
