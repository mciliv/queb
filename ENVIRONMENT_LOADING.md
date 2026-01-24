# Environment Variable Loading System

## Overview

This document explains how environment variables are loaded in the Queb application, with a focus on **loading from the shell environment first**.

## Key Principles

### 🔒 **Shell Environment Takes Priority**
- Environment variables set in the shell/system environment **always take precedence** over `.env` file variables
- This prevents `.env` files from accidentally overriding production settings
- Shell environment variables cannot be overridden by `.env` loading

### 🌐 **Cloud/Production vs Local Development**
- **Production/Cloud**: Uses shell environment variables only (secure, no `.env` loading)
- **Local Development**: Loads `.env` file as fallback for missing variables

### ⚠️ **Critical Variable Validation**
- Validates presence of essential variables (`OPENAI_API_KEY`, `OPENAI_MODEL`)
- Provides clear warnings for missing critical variables
- Helps catch configuration issues early

## Implementation Details

### Environment Detection
The system detects cloud/production environments by checking for:
- `GAE_APPLICATION` (Google App Engine)
- `GOOGLE_CLOUD_PROJECT` (GCP)
- `K_SERVICE` (Cloud Run)
- `FUNCTION_NAME` (Cloud Functions)
- `AWS_EXECUTION_ENV` (AWS Lambda)
- `VERCEL` (Vercel)
- `NETLIFY` (Netlify)
- `NODE_ENV === 'production'`

### Loading Strategy
```javascript
// 1. Detect environment
const isCloud = detectCloudEnvironment();

if (isCloud) {
  // Production: Shell environment only
  console.log('🏭 Using SHELL ENVIRONMENT variables only');
} else {
  // Local dev: Load .env with shell priority
  console.log('🏠 SHELL ENVIRONMENT takes priority, .env provides fallbacks');
  dotenv.config({ path: '.env', override: false }); // override: false = shell wins
}
```

### Priority Order
1. **Shell/System Environment Variables** (highest priority)
2. **`.env` File Variables** (fallback for local development only)
3. **Default Values** (lowest priority)

## Testing

Run the environment loading tests:
```bash
node test-env-loading.js
```

This verifies:
- Shell environment variables take priority
- Cloud detection works correctly
- Critical variable validation functions
- No accidental overrides of shell variables

## Usage Examples

### Local Development
```bash
# Set shell variable (takes priority)
export OPENAI_API_KEY="sk-shell-key"

# .env file has different key (ignored)
# OPENAI_API_KEY=sk-dotenv-key

# Result: OPENAI_API_KEY = "sk-shell-key" (from shell)
```

### Production Deployment
```bash
# Only shell environment variables are used
export OPENAI_API_KEY="sk-production-key"
export NODE_ENV=production

# .env file is completely ignored
```

## Security Benefits

1. **No .env in Production**: Prevents `.env` files from being committed or deployed
2. **Shell Override Protection**: Shell variables can't be accidentally overridden
3. **Clear Separation**: Local vs production configurations are clearly separated
4. **Validation**: Missing critical variables are caught and reported

## Migration Notes

If you're migrating from an old system where `.env` overrode shell variables:
1. Check for any scripts that relied on `.env` overriding shell vars
2. Update those scripts to use shell variables instead
3. The new system is more secure and predictable

---

**Bottom Line**: Environment variables are loaded from the **shell environment first**, with `.env` providing fallbacks only for local development. This ensures production deployments are secure and predictable.