# Molecular Analysis Testing

## Quick Start

```bash
# Run all tests (recommended)
npm test

# Quick molecular test with specific compounds
npm run test:quick water ethanol coffee

# Watch mode for development
npm run test:watch
```

## What `npm test` Does

The unified test runner automatically:

1. **🧪 Unit Tests** - Core functionality validation
2. **🖥️ Visual Interface Tests** - UI components and interactions  
3. **🧬 Molecular Visualization Tests** - 3D molecule rendering
4. **🔗 API Integration Tests** - AI analysis pipeline (if API key available)

## Key Features

- **Seamless Tab Management** - Automatically reuses existing Chrome tabs
- **No Manual Setup** - Everything runs automatically
- **Visual Screenshots** - Captures test results for inspection
- **Parallel Execution** - Runs tests efficiently
- **Smart Fallbacks** - Gracefully handles missing dependencies

## Test Results

Screenshots are saved to `test/screenshots/`:
- Interface screenshots show UI state
- Molecular screenshots show 3D visualizations
- All tests leave Chrome open for manual inspection

## Development Workflow

```bash
# Start development (runs tests first)
npm run dev

# Just run tests
npm test

# Quick test specific molecules
npm run test:quick "red wine" salt caffeine
```

The testing system integrates seamlessly with your development workflow - no complex setup required!