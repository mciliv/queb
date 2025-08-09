# 📸 Screenshot System for LLM Analysis

## Overview
Automated Chrome screenshot capture system that creates images accessible to LLMs for visual analysis of your molecular app.

## Components Created

### 1. Screenshot Service (`backend/services/screenshot-service.js`)
- **Purpose**: Automated Chrome/Puppeteer screenshot capture
- **Features**:
  - Capture basic app state
  - Capture app with text input
  - Capture after analysis is triggered
  - Automatic browser management
  - Screenshot cleanup utilities

### 2. API Endpoints (`backend/api/server.js`)
- **`POST /api/capture-screenshot`** - Basic app screenshot
- **`POST /api/capture-with-input`** - Screenshot with text input
- **`POST /api/capture-analysis`** - Screenshot after triggering analysis
- **`GET /api/screenshots`** - List all screenshots
- **`GET /api/screenshot/:filename`** - Serve screenshot files
- **`GET /api/screenshot-info/:filename`** - Get screenshot metadata
- **`POST /api/cleanup-screenshots`** - Clean up old screenshots

### 3. CLI Tool (`scripts/capture-screenshot.js`)
Easy command-line interface for screenshot capture:
```bash
# Basic app screenshot
node scripts/capture-screenshot.js

# Screenshot with input text
node scripts/capture-screenshot.js --input "water"

# Screenshot after analysis
node scripts/capture-screenshot.js --analysis "caffeine"

# List all screenshots
node scripts/capture-screenshot.js --list
```

## Usage Examples

### Via API (for automated systems)
```javascript
// Capture basic screenshot
const response = await fetch('http://localhost:3000/api/capture-screenshot', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({ filename: 'my-screenshot.png' })
});

// Capture with input
const response = await fetch('http://localhost:3000/api/capture-with-input', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({ 
    inputText: 'aspirin',
    filename: 'aspirin-input.png' 
  })
});
```

### Via CLI (for manual testing)
```bash
# Quick screenshot of current app state
node scripts/capture-screenshot.js

# Test specific molecule input
node scripts/capture-screenshot.js --input "ethanol"

# Capture full analysis workflow
node scripts/capture-screenshot.js --analysis "vitamin C"
```

## LLM Integration

Screenshots are automatically available for LLM analysis via HTTP URLs:

1. **Storage Location**: `/screenshots/` directory
2. **Web Access**: `http://localhost:3000/api/screenshot/[filename]`
3. **Metadata API**: `http://localhost:3000/api/screenshot-info/[filename]`

### Example LLM Workflow
```javascript
// 1. Capture screenshot
const result = await captureScreenshot('water');

// 2. LLM can now access via URL
const screenshotUrl = `http://localhost:3000${result.screenshot.url}`;

// 3. Pass URL to LLM for visual analysis
const analysis = await llm.analyzeImage(screenshotUrl);
```

## Key Features

✅ **Automated Browser Management** - No manual Chrome setup needed
✅ **Multiple Capture Modes** - Basic, input, and analysis screenshots  
✅ **HTTP Accessible** - Screenshots served via web URLs for LLM access
✅ **Metadata Support** - File info, timestamps, and context available
✅ **Cleanup Utilities** - Automatic old screenshot removal
✅ **CLI Interface** - Easy manual testing and automation
✅ **Error Handling** - Robust error recovery and reporting

## File Structure
```
screenshots/           # Screenshot storage (auto-created)
├── app-[timestamp].png
├── input-[name]-[timestamp].png
└── analysis-[name]-[timestamp].png

backend/services/
└── screenshot-service.js   # Core screenshot functionality

scripts/
└── capture-screenshot.js   # CLI utility

backend/api/server.js       # API endpoints (added)
```

## Usage Scenarios

### 1. Visual Testing
Capture app states for automated visual regression testing.

### 2. LLM Analysis
Enable LLMs to "see" your molecular visualizations and provide feedback.

### 3. Documentation
Automatically generate screenshots for documentation and tutorials.

### 4. Debugging
Capture problematic states for issue diagnosis.

## Next Steps

1. **Start the development server**: `./dev`
2. **Test screenshot capture**: `node scripts/capture-screenshot.js --help`
3. **Integrate with LLM**: Use the HTTP URLs for visual analysis
4. **Automate workflows**: Build on the API endpoints for custom automation

The system is now ready for LLM visual analysis integration! 🚀

