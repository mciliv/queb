# Visual Tests

This directory contains automated visual tests for molecular visualization.

## Structure

- **`visual-test-runner.js`** - Main test runner that executes preset visual tests
- **`README.md`** - This file

## Running Tests

### Standalone Runner

Run the visual test runner directly:

```bash
node tests/suites/visual/visual-test-runner.js
```

### Via Jest Integration Tests

Run the integration tests that use the visual test system:

```bash
npm test tests/suites/integration/visual-test-integration.test.js
```

### Unit Tests for Helper Functions

Run unit tests for the underlying helper functions:

```bash
npm test tests/suites/unit/visual-test-helpers.test.js
```

## Test Data

Test data is defined in `src/client/core/constants.js`:

- **`PRESET_VISUAL_TESTS`** - Array of test configurations with labels and SMILES lists
- **`SMILES_NAME_MAP`** - Mapping of SMILES notation to human-readable names

## Architecture

The visual test system is split into:

1. **Pure functions** (`src/client/utils/visual-test-helpers.js`)
   - Path sanitization
   - SDF path resolution
   - Viewer object creation
   - Fully unit-tested

2. **Test runner** (`tests/suites/visual/visual-test-runner.js`)
   - Orchestrates test execution
   - Handles SDF generation
   - Validates results
   - Provides CLI output

3. **Integration tests** (`tests/suites/integration/visual-test-integration.test.js`)
   - End-to-end testing
   - Server interaction validation
   - Error handling verification

## Adding New Tests

Edit `src/client/core/constants.js` and add to `PRESET_VISUAL_TESTS`:

```javascript
{
  label: 'My Test',
  smilesList: ['CCO', 'O', 'C=O']
}
```

Add corresponding name mappings to `SMILES_NAME_MAP` if desired.

## Design Principles

- **Production code is test-free**: No test logic in `App.jsx` or other production components
- **Testable utilities**: All logic extracted into pure, testable functions
- **Multiple test layers**: Unit tests for helpers, integration tests for the full flow
- **Clear separation**: Tests live in `tests/`, utilities in `src/`

