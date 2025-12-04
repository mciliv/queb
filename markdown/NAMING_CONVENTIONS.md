# Naming Conventions

## Files

| Type | Convention | Example |
|------|------------|---------|
| Services/Utils | kebab-case | `molecular-processor.js` |
| React Components | PascalCase | `MoleculeViewer.jsx` |
| Python | snake_case | `sdf_generator.py` |
| Tests | `[name].test.js` | `molecular-processor.test.js` |

## Code

```javascript
// Constants (UPPER_SNAKE_CASE)
const MAX_MOLECULES = 100;

// Variables (camelCase, boolean prefix: is/has/should)
const moleculeData = {};
const isLoading = false;

// Functions (verb+noun, boolean: is/has/can prefix)
function analyzeMolecule() {}
function isValidSmiles() {}
function handleClick() {}

// React components (PascalCase)
function MoleculeViewer() {}
```

## API Endpoints

```
POST   /api/molecules           # Create
GET    /api/molecules/:id       # Read
POST   /api/analyze-text        # Action (verb-noun)
```

## Database

```sql
-- Tables: snake_case, plural
CREATE TABLE analysis_results (...);

-- Columns: snake_case
molecule_name, created_at, user_id
```

## CSS (BEM)

```css
.molecule-viewer { }
.molecule-viewer__canvas { }
.molecule-viewer--loading { }
```

## Git Branches

```
feature/add-protein-support
fix/memory-leak-viewer
refactor/cleanup-services
```

## Guidelines

1. Be descriptive: `processChemicalData()` not `process()`
2. Avoid abbreviations: `temperature` not `temp`
3. Use domain terms: `smiles`, `molecule`, `compound`
4. Indicate purpose: `validateSmiles()` not `check()`
