# Development Guide

## Quick Start

### Chrome Auto-Reload Development (Recommended)
```bash
npm run dev
```
This will:
- Start the application server
- Automatically open Chrome to the app
- Watch for file changes and reload intelligently
- Handle frontend and backend changes appropriately

### Alternative Development Commands

```bash
# Server only (no Chrome integration)
npm run dev:server

# Server with file watching (no auto-reload)
npm run dev:watch

# Chrome auto-reload (explicit)
npm run dev:chrome

# Build frontend and watch for changes
npm run build:watch

# Full development with concurrent frontend building
npm run dev:full
```

## Chrome Auto-Reload System

The new development system replaces live-reload with Chrome DevTools Protocol integration:

### ✅ Benefits
- **Automatic Chrome launch** - No manual browser opening
- **Smart reloading** - Different reload strategies for different file types
- **Better performance** - Uses native Chrome protocols
- **Robust error handling** - Graceful recovery from issues

### 🔄 Reload Strategies
- **Frontend changes** → Soft reload (fastest)
- **Backend changes** → Page reload (medium)
- **Config changes** → Hard reload with cache clear (thorough)

### 📁 Watched Directories
- `src/client/` - Frontend React components
- `src/core/` - Core modules and utilities  
- `src/server/` - Backend API and services

## Application URLs

- **HTTP**: http://localhost:8080
- **HTTPS**: https://localhost:3001
- **Health Check**: http://localhost:8080/health

## File Structure

```
src/
├── client/          # Frontend React application
├── core/            # Shared core modules (Configuration, ErrorHandler, etc.)
├── server/          # Backend API and services
└── ...

scripts/
├── dev-chrome.js    # Chrome auto-reload development script
└── ...

docs/
├── CHROME_AUTO_RELOAD.md    # Detailed Chrome auto-reload documentation
└── ...
```

## Core Modules (Philosophy of Software Design)

The project follows "A Philosophy of Software Design" principles with deep modules:

### Configuration.js
```javascript
const config = require('./src/core/Configuration');
const port = config.get('port');
const dbConfig = config.getDatabaseConfig();
```

### ErrorHandler.js  
```javascript
const errorHandler = require('./src/core/ErrorHandler');
const result = errorHandler.handle(error, context);
```

### PromptEngine.js
```javascript
const promptEngine = require('./src/core/PromptEngine');
const prompt = promptEngine.generateChemicalPrompt(description);
```

### MolecularAnalysisService.js
```javascript
const analysisService = require('./src/core/MolecularAnalysisService');
const result = await analysisService.analyzeText(description);
```

## Testing

```bash
# Run all tests
npm test

# Run specific test suite
npm run test:quick

# Run with coverage
npm test -- --coverage
```

## Building

```bash
# Build frontend once
npm run build

# Build and watch for changes
npm run build:watch
```

## Environment Variables

Create a `.env` file in the project root:

```env
# Required
OPENAI_API_KEY=your_openai_key_here

# Optional
OPENAI_MODEL=gpt-4o
PORT=8080
HTTPS_PORT=3001

# Database (optional)
DB_ENABLED=false
DB_HOST=localhost
DB_PORT=5432
DB_NAME=mol_users
DB_USER=mol_user
DB_PASSWORD=mol_password

# Payments (optional)
PAYMENTS_ENABLED=false
PAYMENTS_DEV_MODE=true
```

## Troubleshooting

### Chrome Not Opening
1. Install Google Chrome
2. Check port 9222 availability
3. Try: `pkill -f "Google Chrome"`

### Server Won't Start
1. Check port 8080 availability
2. Verify environment variables
3. Try: `npm run clean`

### Build Issues
1. Clear build cache: `rm -rf src/client/dist`
2. Reinstall dependencies: `rm -rf node_modules && npm install`
3. Check Node.js version compatibility

### File Changes Not Detected
1. Check watch directory permissions
2. Verify file patterns in DevServer config
3. Restart development server

## Migration from Live-Reload

The Chrome auto-reload system replaces the previous live-reload implementation:

### ✅ Automatic Migration
- No code changes required
- Same npm scripts work
- Better performance and reliability

### ⚠️ Requirements
- Google Chrome must be installed
- Port 9222 must be available for Chrome DevTools

### 🔄 Differences
- Chrome-specific (no Firefox/Safari auto-reload)
- More intelligent reload strategies
- Better integration with development tools

For detailed information about the Chrome auto-reload system, see [docs/CHROME_AUTO_RELOAD.md](docs/CHROME_AUTO_RELOAD.md).

## Additional Documentation

- [Chrome Auto-Reload System](docs/CHROME_AUTO_RELOAD.md) - Detailed guide to the new development system
- [Philosophy of Software Design Improvements](PHILOSOPHY_OF_SOFTWARE_DESIGN_IMPROVEMENTS.md) - Architecture improvements applied
