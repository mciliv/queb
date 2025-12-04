# Queb Hook Architecture

## Overview

Queb's client hooks are designed around **user intent** rather than technical implementation. Each hook represents a specific action a user wants to perform, making the code self-documenting and easy to understand.

## Design Principles

### 1. User-Focused Naming
Hooks are named after what users want to do, not how it's implemented:
- ✅ `useAnalyzeText` - "I want to analyze text"
- ❌ `useApi` - Technical, unclear intent

### 2. Single Responsibility
Each hook does one thing well:
- `useAnalyzeText` - Text analysis only
- `useAnalyzeCamera` - Camera analysis only
- `useMolecule3D` - 3D structure generation only

### 3. Built-in State Management
Hooks manage their own state, exposing only what's needed:
```javascript
const { analyze, molecules, isAnalyzing, error } = useAnalyzeText();
// No need to manage loading states manually
```

## Available Hooks

### Analysis Hooks

#### `useAnalyzeText`
Analyze text descriptions to find molecular composition.

```javascript
const { 
  analyze,      // Function to analyze text
  molecules,    // Result: array of molecules found
  isAnalyzing,  // Loading state
  error,        // Error message if any
  clear         // Clear results
} = useAnalyzeText();

// Usage
const result = await analyze('coffee', { mode: 'database' });
```

#### `useAnalyzeImage`
Analyze uploaded images or photos.

```javascript
const { 
  analyzeImage,     // Analyze image file
  analyzeImageUrl,  // Analyze image from URL
  molecules,        // Found molecules
  identifiedObject, // What was identified
  isAnalyzing,
  error
} = useAnalyzeImage();

// Usage
const result = await analyzeImage(fileInput.files[0]);
```

#### `useAnalyzeCamera`
Live camera analysis with click-to-analyze.

```javascript
const {
  startCamera,   // Start camera stream
  stopCamera,    // Stop camera stream
  analyzeClick,  // Analyze at click coordinates
  cameraActive,  // Is camera running?
  molecules,
  isAnalyzing
} = useAnalyzeCamera();

// Usage
await startCamera(videoElement);
const result = await analyzeClick(x, y);
```

### Data Hooks

#### `useChemicalData`
Look up detailed information about specific chemicals.

```javascript
const {
  lookupChemical,  // Look up by name
  lookupByCAS,     // Look up by CAS number
  lookupBySMILES,  // Look up by SMILES notation
  chemicalInfo,    // Detailed chemical data
  isLoading
} = useChemicalData();

// Usage
const info = await lookupChemical('caffeine');
// Returns: formula, IUPAC name, properties, etc.
```

#### `useMolecule3D`
Generate and manage 3D molecular structures.

```javascript
const {
  generate3DFromName,  // Generate from chemical name
  generateMultiple3D,  // Generate multiple at once
  getViewerData,      // Get data for 3D viewer
  structures,         // Generated structures
  isGenerating
} = useMolecule3D();

// Usage
const structure = await generate3DFromName('aspirin');
```

## Migration Guide

### Before (Generic useApi)
```javascript
const { structuresFromText, generateSDFs } = useApi();

const handleAnalyze = async (text) => {
  setLoading(true);
  try {
    const result = await structuresFromText(text, 'foodb');
    const molecules = result.molecules || [];
    
    if (molecules.length > 0) {
      const smiles = molecules.map(m => m.smiles);
      const sdfs = await generateSDFs(smiles);
      // Complex mapping...
    }
    setResults(molecules);
  } catch (err) {
    setError(err.message);
  } finally {
    setLoading(false);
  }
};
```

### After (Focused Hooks)
```javascript
const { analyze, molecules, isAnalyzing, error } = useAnalyzeText();
const { generateMultiple3D } = useMolecule3D();

const handleAnalyze = async (text) => {
  const result = await analyze(text, { mode: 'database' });
  if (result.molecules.length > 0) {
    await generateMultiple3D(result.molecules);
  }
};
```

## Benefits

1. **Clearer Code**: Function names express intent
2. **Less Boilerplate**: State management is built-in
3. **Better Errors**: Consistent error handling
4. **Easier Testing**: Each hook can be tested independently
5. **Better Types**: Each hook has specific TypeScript types

## Implementation Details

### Core API Hook
All user-focused hooks use `useApiCore` internally for HTTP communication:

```javascript
// useApiCore handles:
// - HTTP requests
// - Retries
// - Timeouts
// - Cancellation

// User hooks handle:
// - State management
// - Data transformation
// - User-friendly APIs
```

### Caching Strategy
Some hooks implement caching to improve performance:

```javascript
// useMolecule3D caches generated structures
const structures = new Map(); // name -> structure data

// useAnalyzeText could cache recent analyses
const recentAnalyses = new Map(); // text -> results
```

## Best Practices

1. **Use the Right Hook**: Choose based on user action, not data type
2. **Let Hooks Manage State**: Don't duplicate state management
3. **Handle Errors Gracefully**: All hooks expose error states
4. **Compose When Needed**: Combine hooks for complex workflows

```javascript
// Good: Composing hooks for a complete flow
function AnalyzeAndVisualize() {
  const { analyze, molecules } = useAnalyzeText();
  const { generateMultiple3D, structures } = useMolecule3D();
  
  const handleFullAnalysis = async (text) => {
    const result = await analyze(text);
    if (result.molecules.length > 0) {
      const structures = await generateMultiple3D(result.molecules);
      // Now you have both analysis and 3D structures
    }
  };
}
```

## Future Hooks

Potential hooks to add based on user needs:

- `useAnalyzeVideo` - Analyze video streams
- `useMoleculeSearch` - Search molecular databases
- `useMoleculeComparison` - Compare molecules
- `useAnalysisHistory` - Track analysis history
- `useMoleculeExport` - Export molecules in various formats


