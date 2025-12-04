# Contributing

## Code Structure

### Frontend (`src/client/`)
- `components/` - React components (App.jsx, MoleculeViewer.jsx)
- `hooks/` - Custom hooks (useApi.js)
- `utils/` - Frontend utilities
- `assets/` - Static assets

### Backend (`src/server/`)
- `api/` - Express endpoints (server.js)
- `services/` - Business logic (Structuralizer.js, molecular-processor.js)
- `utils/` - Backend utilities

### Core (`src/core/`)
- `PromptEngine.js` - AI prompt templates
- `ErrorHandler.js` - Error handling
- `ServiceContainer.js` - Dependency injection

## Patterns

### Adding a New Input Mode
```javascript
// src/client/components/NewInputMode.jsx
export default function NewInputMode({ onAnalyze }) {
  const handleSubmit = (data) => {
    onAnalyze({ type: 'new-mode', data });
  };
  return <div className="input-mode">{/* UI */}</div>;
}
```

### Adding an API Endpoint
```javascript
// In src/server/api/server.js
router.post('/new-endpoint', async (req, res, next) => {
  try {
    const { data } = req.body;
    if (!data) return res.status(400).json({ error: 'Data required' });
    const result = await processData(data);
    res.json({ success: true, result });
  } catch (error) {
    next(error);
  }
});
```

## Naming Conventions

| Type | Convention | Example |
|------|------------|---------|
| Files | kebab-case | `molecular-processor.js` |
| Components | PascalCase | `MoleculeViewer.jsx` |
| Functions | camelCase | `processChemical()` |
| Constants | UPPER_CASE | `MAX_MOLECULES` |

## Commit Messages

```
feat: add support for protein structures
fix: resolve memory leak in molecule viewer
docs: update API documentation
refactor: simplify chemical name resolution
test: add unit tests for name resolver
```

## Best Practices

**Performance:**
- Cache expensive operations
- Minimize API calls
- Use lazy loading

**Security:**
- Never commit secrets
- Validate all inputs
- Sanitize user content

**Testing:**
- 80%+ code coverage
- Test edge cases
- Mock external APIs
