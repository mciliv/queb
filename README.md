# Estimate & visualize what molecules are contained in a given section(s) of space

Upload an image, take a photo, or describe an object → Get molecules → See interactive 3D molecular structures

**Example**: Photo of coffee → Identifies caffeine, water, etc. → Displays 3D molecular models

UI structure:
# Text input (⌘k)
# camera|image|link input buttons in which 0/1 can be selected
## Camera 
### Mobile
    Rectangle reticle in center of screen, click screen to capture
### Non-mobile
    Click area in the image in which the click position is the center of a red square, whose crop is used as input for identification

# Molecule x object grid
upon entered input, 
Options:
1/Single colmn
2/a column is added to the right if a column already exists
A column contains a header, specified by the object specification or inference (Load immediately upon input submittance)
each molecule have name title which links to wiki 
3dmoljs gridviewer seems best

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