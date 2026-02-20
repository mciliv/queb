# Queb — Chemical Structure Analyzer

Point a camera at something, upload a photo, paste a link, or type a name → AI identifies the object and its chemical compounds → Interactive 3D molecular structures appear.

**Example**: Photo of coffee → Caffeine, chlorogenic acid, trigonelline → Rotatable 3D models

## Quick Start

```bash
npm install
cp .env.example .env   # add your API key
npm start              # http://localhost:8080
```

## Environment Variables

```bash
# Required — choose one provider
AI_PROVIDER=openai
OPENAI_API_KEY=sk-...

# or
AI_PROVIDER=xai
XAI_API_KEY=xai-...

# Optional
PORT=8080
DATABASE_ENABLED=false   # set true + DATABASE_URL to enable user auth
```

## Input Methods

| Mode | How |
|------|-----|
| Text | Type an object name ("coffee", "aspirin") |
| Camera | Live feed — click anything to analyze it |
| Photo | Upload an image |
| Link | Paste an image URL or drag one in |

## Architecture

```
src/
├── server/
│   ├── api/
│   │   └── server.js          # Entry point (npm start)
│   │   └── app.js             # Express app + route wiring
│   │   └── routes/
│   │       ├── chemical-analysis.js
│   │       └── hotel.js
│   └── services/
│       ├── AIService.js        # @ai-sdk unified wrapper (OpenAI / xAI)
│       ├── molecular-processor.js  # SDF generation (PubChem → Python/RDKit)
│       ├── name-resolver.js    # Name → CID → SMILES (PubChem, ChEMBL, CACTUS)
│       ├── structuralizer.js   # Core prediction pipeline
│       └── database-recommender.js
└── client/
    ├── components/
    │   ├── App.jsx             # Root component
    │   ├── CameraSection.jsx
    │   ├── PhotoSection.jsx
    │   ├── LinkSection.jsx
    │   ├── TextInput.jsx
    │   ├── MoleculeViewer.jsx  # 3Dmol.js wrapper
    │   └── MolecularColumn.jsx
    └── hooks/
        └── useApi.js

core/
├── services.js        # DI container factory
├── ServiceContainer.js
└── PromptEngine.js
```

**Entry point**: `src/server/api/server.js`
**DI container**: `src/core/services.js` — ~15 lazily-resolved services
**AI layer**: `AIService.js` uses Vercel AI SDK (`@ai-sdk/openai`, `@ai-sdk/xai`) — provider is switched via `AI_PROVIDER` env var with no code changes
**SDF generation**: PubChem download first; falls back to local Python/RDKit subprocess
**Database**: Optional PostgreSQL (`pg`) for user accounts; app runs fully without it

## API

### POST /api/structuralize
```json
{ "text": "coffee", "lookupMode": "GPT-5" }
→ { "object": "coffee", "chemicals": [{ "name": "Caffeine", "smiles": "...", "sdfPath": "..." }] }
```

### POST /api/structuralize-image
```json
{ "imageBase64": "<base64>", "x": 150, "y": 200 }
```

### POST /api/generate-sdfs
```json
{ "smiles": ["CN1C=NC2=C1C(=O)N(C(=O)N2C)C"] }
→ { "sdfPaths": ["/sdf_files/caffeine.sdf"] }
```

### POST /api/find-chemical-contents
```json
{ "item": "green tea" }
→ { "item": "green tea", "selectedDatabase": "...", "chemicals": [...] }
```

### GET /api/3d-structure/:name
Returns SDF content + metadata from PubChem for direct 3D display.

## Development

```bash
npm run dev          # server + frontend watch (requires concurrently)
npm test             # unit tests
npm run test:integration
```

## Tech Stack

- **Frontend**: React 19, 3Dmol.js
- **Backend**: Node.js, Express
- **AI**: Vercel AI SDK — OpenAI GPT-4o Vision or xAI Grok
- **Chemistry**: PubChem REST API, ChEMBL, CACTUS, RDKit (optional)
- **3D**: SDF files served from `/sdf_files/`, rendered client-side

## License

MIT
