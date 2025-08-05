# Background Testing System

Automated testing system that intelligently runs tests based on breaking change probability with persistent, reliable execution across machine restarts.

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

## Persistent Testing Setup

### Installation Methods

#### 1. PM2 (Process Manager)
```bash
# Install PM2 globally
npm install -g pm2

# Start test runners with PM2
npm run pm2:start

# Enable startup on boot
pm2 startup
pm2 save

# Check status
npm run pm2:status
```

**Features:**
- Automatic restart on crash
- Memory limit protection
- Log management
- Process monitoring
- Scheduled restarts to prevent memory leaks

#### 2. Systemd Service (Linux)
```bash
# Run installation script as root
sudo bash scripts/install-test-runner.sh

# Choose option 2 for systemd
# Service will start automatically on boot

# Check service status
systemctl status mol-test-runner

# View logs
journalctl -u mol-test-runner -f
```

**Features:**
- Native OS integration
- Resource limits (CPU/Memory)
- Security hardening
- Automatic startup on boot
- System journal logging

#### 3. Docker Container
```bash
# Start test runner container
npm run docker:test

# Check logs
npm run docker:test:logs

# Stop container
npm run docker:test:stop
```

**Features:**
- Isolated environment
- Consistent across platforms
- Easy deployment
- Built-in health checks
- Resource constraints

### Scheduled Testing

The test scheduler runs automatically with:
- **Every 30 minutes**: Smoke tests
- **Every 2 hours**: Backend unit tests
- **Daily at 3 AM**: Full test suite
- **Every 5 minutes**: Health checks

### Health Monitoring

Health checks monitor:
- Disk space (< 90% usage)
- Memory usage (< 85% heap)
- Git status (< 50 uncommitted files)
- Node modules availability
- Test failure rate (< 20%)

### Log Management

Logs are automatically rotated:
- Daily rotation
- 7 days retention
- Compression after 1 day
- Located in `/var/log/mol/` or `./logs/`

### Reliability Features

1. **Multiple Layers**:
   - PM2 process management
   - Systemd service supervision
   - Cron job fallback
   - Docker container isolation

2. **Auto-Recovery**:
   - Automatic restart on crash
   - Memory limit restart
   - Scheduled process refresh
   - Health check recovery

3. **Monitoring**:
   - Real-time process status
   - Test execution statistics
   - System health metrics
   - Failure notifications

### Quick Start Guide

```bash
# Option 1: Development setup with PM2
npm install
npm run pm2:start
pm2 startup
pm2 save

# Option 2: Production setup with systemd
sudo bash scripts/install-test-runner.sh
# Choose option 2

# Option 3: Container setup
docker-compose -f docker-compose.test.yml up -d

# Check everything is running
npm run test:status
npm run pm2:status
```

### Troubleshooting Persistent Tests

**Tests not starting after reboot:**
```bash
# Check PM2 startup
pm2 startup
pm2 save

# Check systemd
systemctl enable mol-test-runner
systemctl start mol-test-runner

# Check Docker
docker ps -a
docker-compose -f docker-compose.test.yml up -d
```

**High memory usage:**
```bash
# PM2 memory management
pm2 restart mol-test-scheduler

# Check memory limits
pm2 show mol-test-scheduler | grep memory
```

**Missing test results:**
```bash
# Check health status
node scripts/health-check.js

# View recent logs
tail -n 100 logs/test-watcher-out.log

# Check scheduler status
pm2 logs mol-test-scheduler --lines 50
```