# Unified Testing System

## Quick Start

Run all tests:
```bash
npm test
```

Run specific test types:
```bash
npm run test:quick         # Quick molecular tests
npm run test:visual        # Visual interface tests  
npm run test:molecular     # Molecular visualization tests
npm run test:pipeline      # Full pipeline tests (requires API key)
npm run test:persistence   # Tab management tests
npm run test:demo          # Visual demo tests (slow/visible)
npm run test:cli           # Command line molecular tests
npm run test:fake-data     # Fake data injection tests
```

## How It Works

The unified test runner (`test/scripts/unified-test-runner.js`) orchestrates all test types:

### Core Test Types
1. **Unit Tests** - Core functionality validation
2. **Visual Tests** - UI interaction and display  
3. **Molecular Tests** - End-to-end molecular analysis
4. **Pipeline Tests** - Full AI → SDF → Visualization flow
5. **Persistence Tests** - Chrome tab management
6. **Injection Tests** - Visual data injection

### Utility Tests
- **Demo Tests** - Slow, visible test demonstrations
- **CLI Tests** - Command line molecular testing
- **Fake Data Tests** - Backend pipeline without AI

## Features

- **Tab Persistence** - Reuses Chrome tabs across test runs
- **Smart Skipping** - Skips API tests when no key available  
- **Parallel Execution** - Runs compatible tests simultaneously
- **Screenshot Capture** - Visual test results saved to `test/screenshots/`
- **Comprehensive Reporting** - Clear pass/fail status for each test type
- **Organized Structure** - All test files properly organized in `test/` directory

## Test File Organization

```
test/
├── scripts/           # Standalone test scripts
│   ├── unified-test-runner.js    # Main test orchestrator
│   ├── test-molecular-ui.js      # CLI molecular tests
│   ├── test-visual-demo.js       # Visual demo tests
│   └── test-inject-fake-data.js  # Fake data tests
├── integration/       # Jest integration tests
│   ├── visual-interface.test.js
│   ├── persistent-tab-tests.test.js
│   ├── full-pipeline-visualization.test.js
│   └── molecular-accuracy.test.js
└── utils/            # Test utilities
    ├── auto-tab-connector.js
    └── chrome-tab-manager.js
```

## Arguments

- No args: Full comprehensive test suite
- `water ethanol`: Quick tests with specific molecules
- `--quick [molecules]`: Fast molecular visualization tests
- `--visual`: Visual interface and injection tests
- `--molecular`: Molecular visualization tests only
- `--pipeline`: Full pipeline tests (requires API key)
- `--persistence`: Tab management tests
- `--demo`: Slow/visible demo tests
- `--cli`: Command line tests
- `--fake-data`: Backend pipeline tests without AI

## Development Workflow

```bash
# Start development (runs tests first)
npm run dev

# Just run tests
npm test

# Quick test specific molecules
npm run test:quick water ethanol coffee

# Run visual tests only
npm run test:visual

# Run demo (visible, slow tests)
npm run test:demo
```

The testing system integrates seamlessly with your development workflow - everything is organized under one unified runner!