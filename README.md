# Molecular Analysis App

AI-powered molecular identification and 3D visualization from text, images, or camera input.

## Features

- **Text Input**: Describe molecules or chemical compounds
- **Camera/Image**: Upload photos or take pictures for molecular analysis
- **3D Visualization**: Interactive molecular models using 3Dmol.js
- **AI-Powered**: OpenAI Vision API for chemical identification

## Quick Start

### Development
```bash
# Generic development (no domain)
../web/dev serve

# With queb profile
../web/dev serve
# or
PROFILE=queb ../web/dev serve
```

### Deployment
```bash
# Generic deployment
../web/dev deploy

# Queb-specific deployment
../web/dev deploy
# or
PROFILE=queb ../web/dev deploy
```

### Status Check
```bash
# Generic status
../web/dev status

# Queb-specific status
../web/dev status
# or
PROFILE=queb ../web/dev status
```

## Profiles

The app supports domain-specific profiles:

- **Generic**: No domain configuration (localhost only)
- **Queb**: `config/profiles/queb.env` (queb.space domain)

### Using Profiles

1. **Environment Variable**: `PROFILE=queb ../web/dev serve`
2. **Command Line**: `../web/dev --profile=queb serve`
3. **Profile Scripts**: `../web/dev serve`, `../web/dev deploy`, `../web/dev status`

## Configuration

- **Environment**: `.env` file for API keys
- **Profiles**: `config/profiles/[name].env` for domain settings
- **Project**: `config/project.js` for project metadata
- **SSL**: Auto-generated certificates for local development

## SSL Certificates

Generate trusted SSL certificates for local development:

```bash
# Generate certificates (stored in config/certs/)
../web/generate-certs.js

# Generate with custom domains
../web/generate-certs.js myapp.local dev.myapp.local

# Force regeneration
../web/generate-certs.js --force
```

### Prerequisites

Install [mkcert](https://github.com/FiloSottile/mkcert) and run:
```bash
mkcert -install
```

## Requirements

- Node.js
- Python 3.12+
- Google Cloud CLI (for deployment)
- OpenAI API Key (in .env file)