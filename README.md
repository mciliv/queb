# Molecular Space Analyzer

A web application that analyzes the molecular composition of spaces, objects, and images using AI vision and 3D molecular visualization.

## What It Does

Upload an image, take a photo, or describe an object → Get AI-identified molecules → See interactive 3D molecular structures

**Example**: Photo of coffee → Identifies caffeine, water, etc. → Displays 3D molecular models

## Quick Start

```bash
npm run dev
```

Visit `https://localhost:3001` or `http://localhost:8080`

## Project Structure

```
molecular-space-analyzer/
├── frontend/                    # Web application interface
│   ├── core/                   # Core application logic
│   │   ├── molecular-space-app.js  # Main application class
│   │   ├── modular-molecular-app.js # Modular architecture
│   │   └── system-health-checker.js # Health monitoring
│   ├── components/             # UI components
│   │   ├── input/             # Input handling (camera, text, photo)
│   │   ├── visualization/     # 3D molecular display
│   │   └── ui/               # General UI elements
│   └── assets/               # Static assets (icons, styles)
├── backend/                   # Server API
│   ├── api/                  # Express routes & servers
│   ├── services/             # Business logic
│   │   ├── AtomPredictor.js  # AI vision analysis
│   │   └── molecular-processor.js # SMILES/SDF processing
│   └── schemas/              # Data validation
├── molecular-conversion/      # Chemical format conversion
│   └── processors/
│       └── sdf.py           # SMILES → SDF conversion
├── test/                     # Organized test suites
│   ├── suites/              # Test files by type
│   │   ├── unit/           # Unit tests
│   │   ├── integration/    # Integration tests
│   │   └── visual/         # Visual/UI tests
│   └── support/            # Test utilities & config
└── archive/                 # Historical research code
    └── molecular-docking-research/ # Archived docking pipeline
```

## Key Features

- **AI Image Analysis** - OpenAI Vision API identifies molecules in photos
- **Text Analysis** - Describe objects to get molecular composition
- **3D Visualization** - Interactive molecular structures using 3DMol.js
- **Real-time Camera** - Point camera at objects for instant analysis
- **SMILES Processing** - Converts chemical notation to 3D models

## Technology Stack

- **Frontend**: Vanilla JavaScript, 3DMol.js, CSS Grid
- **Backend**: Node.js, Express, OpenAI API
- **Chemistry**: RDKit (Python), SMILES notation, SDF files
- **Testing**: Jest, Puppeteer
- **Development**: Live reload, HTTPS, automated testing

## Development Commands

```bash
npm run dev              # Full development server with live reload
npm run test            # Run comprehensive test suite
npm run test:unit       # Unit tests only
npm run test:integration # Integration tests only
npm run clean           # Clean up processes and temporary files
```

## Architecture Highlights

### Focused Codebase
- **Separated Research Code**: Molecular docking research archived separately
- **Single Purpose**: Web app focused solely on molecular space analysis
- **Clear Dependencies**: Only essential chemistry conversion tools included

### Modular Design
- **Component-Based**: Reusable UI components
- **Service Layer**: Separate business logic from API routes
- **Test Organization**: Tests grouped by purpose and scope

### Performance
- **Lazy Loading**: 3D structures load on demand
- **Caching**: SDF files cached to avoid regeneration
- **Optimized Rendering**: Efficient 3D molecular display

## Contributing

The codebase is organized for clarity and focused development:

1. **Frontend changes** → `frontend/` directory
2. **API changes** → `backend/` directory  
3. **Chemical processing** → `molecular-conversion/` directory
4. **Tests** → `test/suites/` organized by type

## License

MIT License - See LICENSE file for details.