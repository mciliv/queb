# Visual Testing

## Purpose

Visual testing runs chemical analysis on predefined inputs and automatically visualizes results in the browser—separating debugging complexity from visualization verification.

## Test Modes

| Mode | Script | Use Case |
|------|--------|----------|
| Real | `npm run test:visual` | Live AI/API integration testing |
| Mock | `npm run test:visual-mock` | Fast UI testing, offline, predictable |

## Test Cases

| Test | Input | Expected Molecules |
|------|-------|-------------------|
| coffee | "coffee" | caffeine, water, chlorogenic acid |
| saltwater | "salt water" | sodium chloride, water |
| aspirin | "aspirin tablet" | acetylsalicylic acid, cellulose |
| apple | "fresh apple" | fructose, glucose, malic acid |

## How It Works

1. **Analysis**: Loads test config → runs Structuralizer → generates SDF files → saves JSON
2. **Visualization**: Opens browser with `?test=[name]&autoload=true` → displays 3D structures

## Test Data Format

```json
{
  "testName": "coffee",
  "object": "coffee",
  "chemicals": [
    { "name": "caffeine", "sdfPath": "/sdf_files/caffeine.sdf", "status": "ok" }
  ],
  "metadata": { "totalMolecules": 3, "successfulSdfs": 2 }
}
```

## Adding Test Cases

Edit `scripts/run-visual-test.js`:

```javascript
const TEST_CASES = {
  mytest: {
    input: "my test object",
    expectedMolecules: ["molecule1", "molecule2"],
    lookupMode: "GPT-5"
  }
};
```

## Best Practices

- Start with simple tests (saltwater)
- Monitor SDF success rates in metadata
- Always verify 3D structures visually
- Clear browser cache between runs if needed
