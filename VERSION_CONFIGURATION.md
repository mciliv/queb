# Version Configuration System

## Overview
The molecular analysis app now supports programmatic control of React vs Vanilla JS versions through a flexible configuration system. This allows you to control which version loads, whether toggling is allowed, and how the system behaves in different environments.

## Configuration Options

### Core Settings
```javascript
{
  defaultVersion: 'vanilla',    // 'vanilla' | 'react'
  allowToggle: true,            // Enable/disable version switching
  showToggleButton: true,       // Show/hide the toggle button
  persistChoice: true,          // Remember user's choice in localStorage
  performanceMode: false        // Force vanilla JS, hide toggle for production
}
```

### Environment-Specific Overrides
The system automatically detects the environment and applies appropriate settings:

- **Development** (`localhost`, `127.0.0.1`, `dev.*`): Full toggle functionality
- **Production** (other domains): Performance optimized, no toggle
- **Testing** (`?test=true` or `test.*` domains): React default, toggle enabled

## Usage Methods

### 1. URL Parameters
Control behavior via URL parameters:
```
http://localhost:3001/?version=react&toggle=false&showToggle=false
```

**Available Parameters:**
- `version`: `vanilla` | `react`
- `toggle`: `true` | `false`
- `showToggle`: `true` | `false`
- `persist`: `true` | `false`

### 2. Configuration Presets
Apply predefined configurations:
```javascript
// Apply production preset
applyPreset('production');

// Apply React-only mode
applyPreset('reactOnly');

// Apply vanilla-only mode
applyPreset('vanillaOnly');
```

**Available Presets:**
- `development`: Full toggle functionality
- `production`: Performance optimized, no toggle
- `testing`: React default, toggle enabled
- `reactOnly`: React only, no toggle
- `vanillaOnly`: Vanilla JS only, no toggle

### 3. Runtime Configuration
Modify configuration at runtime:
```javascript
// Update configuration
setVersionConfig({
  defaultVersion: 'react',
  allowToggle: false,
  showToggleButton: false
});

// Get current configuration
const config = getVersionConfig();

// Force load specific version
forceLoadVersion('react');
```

### 4. Direct Configuration File
Edit `frontend/core/version-config.js` to change default behavior:
```javascript
export const VERSION_CONFIG = {
  defaultVersion: 'vanilla',
  allowToggle: true,
  showToggleButton: true,
  persistChoice: true,
  performanceMode: false
};
```

## Console Commands

### Configuration Management
```javascript
// Get current configuration
getVersionConfig()

// Update configuration
setVersionConfig({ allowToggle: false })

// Apply preset
applyPreset('production')

// Get active config with environment overrides
getActiveConfig()

// Get final config with URL overrides
getFinalConfig()
```

### Version Control
```javascript
// Switch versions
switchToReact()
switchToVanilla()

// Force load version (bypasses restrictions)
forceLoadVersion('react')

// Get current version
getCurrentVersion()
```

### Debugging
```javascript
// Log version information
logVersionInfo()

// Test version switching
testVersionSwitch()
```

## Environment Detection

### Development Environment
- **Hostnames**: `localhost`, `127.0.0.1`, `dev.*`
- **Behavior**: Full toggle functionality, vanilla JS default
- **Use Case**: Development and testing

### Production Environment
- **Hostnames**: Any domain not matching development patterns
- **Behavior**: Performance optimized, no toggle, vanilla JS only
- **Use Case**: Live production deployment

### Testing Environment
- **URL Parameter**: `?test=true`
- **Hostnames**: `test.*`
- **Behavior**: React default, toggle enabled
- **Use Case**: React-specific testing

## Performance Modes

### Performance Mode Enabled
```javascript
{
  performanceMode: true,
  defaultVersion: 'vanilla',
  allowToggle: false,
  showToggleButton: false,
  persistChoice: false
}
```
- Forces vanilla JS for maximum performance
- Disables all toggle functionality
- No localStorage usage
- Ideal for production deployments

### Development Mode
```javascript
{
  performanceMode: false,
  defaultVersion: 'vanilla',
  allowToggle: true,
  showToggleButton: true,
  persistChoice: true
}
```
- Full toggle functionality
- User choice persistence
- Visual toggle button
- Ideal for development and testing

## Integration Examples

### Production Deployment
```javascript
// In your deployment script
applyPreset('production');
// Result: Vanilla JS only, no toggle, maximum performance
```

### React Development
```javascript
// For React-specific development
applyPreset('reactOnly');
// Result: React only, no toggle, React default
```

### A/B Testing
```javascript
// Randomly assign version for A/B testing
const version = Math.random() > 0.5 ? 'react' : 'vanilla';
setVersionConfig({ defaultVersion: version, allowToggle: false });
```

### Feature Flags
```javascript
// Enable React features for specific users
if (user.hasReactFeatures) {
  setVersionConfig({ defaultVersion: 'react', allowToggle: true });
} else {
  setVersionConfig({ defaultVersion: 'vanilla', allowToggle: false });
}
```

## File Structure
```
frontend/
├── core/
│   ├── version-config.js     # Configuration definitions
│   ├── version-loader.js     # Version switching logic
│   ├── app.js               # Vanilla JS implementation
│   ├── App.jsx              # React implementation
│   └── index.html           # Main entry point
└── assets/
    └── style.css            # Version toggle styles
```

## Best Practices

### 1. Environment-Specific Configuration
- Use environment detection for automatic configuration
- Override with URL parameters for testing
- Use presets for common scenarios

### 2. Performance Optimization
- Enable performance mode in production
- Disable toggle functionality for production
- Use vanilla JS as default for better performance

### 3. Development Workflow
- Use development preset for local development
- Enable toggle for easy testing
- Use URL parameters for quick configuration changes

### 4. Testing Strategy
- Use testing preset for React-specific tests
- Use vanilla-only preset for performance tests
- Use URL parameters for automated testing

## Troubleshooting

### Common Issues
1. **Toggle not working**: Check `allowToggle` configuration
2. **Button not visible**: Check `showToggleButton` configuration
3. **Wrong version loading**: Check `defaultVersion` and environment detection
4. **Configuration not applying**: Ensure configuration is loaded before version loader

### Debug Commands
```javascript
// Check current state
logVersionInfo()

// Verify configuration
console.log(getVersionConfig())

// Test environment detection
console.log(getActiveConfig())
``` 