# Molecular Visualization App - Reorganization Plan

## Current Issues
1. **Unclear Purpose**: Generic name "mol" doesn't convey molecular visualization intent
2. **Browser Profile Conflicts**: Multiple Chrome profiles causing singleton lock errors
3. **Scattered Test Structure**: Integration tests mixed with different concerns
4. **Naming Inconsistency**: Files named generically without clear purpose indication

## New Structure & Naming Convention

### 1. Project Rename
- **Current**: `mol` 
- **New**: `molecular-space-analyzer`
- **Purpose**: Clearly indicates the app analyzes molecular composition of spaces/objects

### 2. Directory Structure
```
molecular-space-analyzer/
├── src/                           # Main source code
│   ├── frontend/                  # Client-side application
│   │   ├── components/           # React components
│   │   │   ├── input/           # Input handling (TextInput, CameraInput, PhotoInput)
│   │   │   ├── visualization/   # 3D molecular display (MolecularViewer, ResultsGrid)
│   │   │   ├── analysis/        # Analysis results (ChemicalSummary, ErrorDisplay)
│   │   │   └── ui/             # Generic UI (Layout, Shortcuts, Theme)
│   │   ├── core/               # Core application logic
│   │   │   ├── molecular-app.js # Main app class (renamed from app.js)
│   │   │   ├── api-client.js   # API communication
│   │   │   └── constants.js    # App constants
│   │   └── assets/             # Static assets
│   ├── backend/                # Server-side API
│   │   ├── api/               # Express routes
│   │   ├── services/          # Business logic
│   │   └── utils/             # Utilities
│   └── chemistry/             # Chemical processing logic
├── tests/                     # All test files (renamed from test/)
│   ├── unit/                 # Unit tests
│   ├── integration/          # Integration tests (organized by feature)
│   │   ├── molecular-analysis/ # Molecular analysis workflow tests
│   │   ├── browser-automation/ # Browser-based tests
│   │   └── api-integration/   # API endpoint tests
│   ├── fixtures/             # Test data and mocks
│   ├── utils/               # Test utilities
│   └── profiles/            # Browser profiles (organized to prevent conflicts)
│       ├── molecular-testing/ # Primary testing profile
│       └── persistence-testing/ # Secondary profile for tab persistence tests
├── docs/                    # Documentation
├── scripts/                 # Build and deployment scripts
└── config/                  # Configuration files
```

### 3. Component Renaming
- `app.js` → `molecular-space-analyzer.js` (main app class)
- `Results.jsx` → `MolecularVisualizationResults.jsx`
- `TextInput.jsx` → `MolecularAnalysisInput.jsx`
- Generic test files → Purpose-specific names

### 4. Test Organization
- **Browser Profile Separation**: Different profiles for different test types
- **Feature-based Grouping**: Tests grouped by molecular analysis features
- **Clear Naming**: Each test file indicates its specific purpose

### 5. Fix Browser Profile Conflicts
- Create separate Chrome profiles for different test suites
- Implement proper cleanup between test runs
- Add profile lock detection and cleanup

## Implementation Steps
1. Rename project and update package.json
2. Reorganize directory structure
3. Fix browser profile conflicts in tests
4. Update component names to reflect molecular analysis purpose
5. Update all imports and references
6. Update documentation to reflect new structure