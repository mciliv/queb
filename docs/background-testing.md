# Background Testing System

Automated testing system that intelligently runs tests based on breaking change probability.

## Features

🔍 **Intelligent Change Detection** - Analyzes git changes to determine breaking change probability
🧪 **Smart Test Selection** - Runs appropriate test suites based on changed files  
⚡ **Non-blocking** - Tests run in background without interrupting development
🎯 **Priority-based** - Critical tests (smoke, backend) run first
📊 **Real-time Feedback** - Clear status reporting and logging

## Automatic Triggers

### Git Hooks
- **Pre-commit**: Analyzes staged changes, runs relevant tests before commit
- **Post-commit**: Starts background tests after successful commit

### File Watching  
- **Development Mode**: Watches for file changes during development
- **Continuous Integration**: Monitors git changes for automated testing

## Usage Commands

```bash
# Manual background test run
npm run test:background

# Check last test results  
npm run test:status

# Start continuous file watching
npm run test:watch-continuous

# Development file watcher with immediate feedback
npm run dev:test-watch
```

## Breaking Change Detection

### High Probability (Full Test Suite)
- API endpoints (`backend/api/`, `backend/services/`)
- Schema changes (`backend/schemas/`)
- Infrastructure (`package.json`, `jest.config.js`, `.husky/`)

### Medium Probability (Targeted Tests)
- Backend code changes (`backend/**/*.js`)
- Frontend code changes (`frontend/**/*.js`)
- Test file modifications

### Low Probability (Smoke Tests Only)
- Documentation changes
- Non-critical asset updates

## Test Execution Priority

1. **Smoke Tests** (always run - blocks commit if failed)
2. **Backend Unit Tests** (critical - blocks if failed)  
3. **Frontend Unit Tests** (non-blocking due to mock issues)
4. **Integration Tests** (when API changes detected)

## Log Files

- **Background tests**: `/tmp/mol-background-tests.log`
- **Pre-commit output**: Console during commit process
- **Development watcher**: Real-time console output

## Integration with Development Workflow

### During Development
```bash
# Start file watcher for immediate feedback
npm run dev:test-watch
```

### Before Committing
- Pre-commit hook automatically analyzes changes
- Runs appropriate tests based on change type
- Blocks commit if critical tests fail

### After Committing  
- Post-commit hook starts background tests
- Full results available via `npm run test:status`

## Configuration

Edit `scripts/background-test-runner.js` to customize:
- Breaking change detection patterns
- Test suite priorities  
- Timeout settings
- Log output format

## Troubleshooting

### Tests Not Running
```bash
# Check if hooks are executable
ls -la .husky/

# Test runner status
npm run test:status

# Manual test run
npm run test:background
```

### Hook Issues
```bash
# Reinstall husky hooks
npx husky install
```

### Background Process Management
```bash
# Kill background tests
pkill -f "background-test-runner"

# Check running processes
ps aux | grep test
```

## Best Practices

1. **Let hooks work** - Don't bypass pre-commit checks
2. **Monitor background results** - Check `npm run test:status` regularly
3. **Use dev watcher** - Run `npm run dev:test-watch` during active development
4. **Fix critical failures** - Backend test failures should be addressed immediately
5. **Mock issues** - Frontend test failures are currently non-blocking due to known mocking issues