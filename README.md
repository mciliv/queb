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
./scripts/dev serve

# With queb profile
./scripts/dev-queb serve
# or
PROFILE=queb ./scripts/dev serve
```

### Deployment
```bash
# Generic deployment
./scripts/dev deploy

# Queb-specific deployment
./scripts/deploy-queb
# or
PROFILE=queb ./scripts/dev deploy
```

### Status Check
```bash
# Generic status
./scripts/dev status

# Queb-specific status
./scripts/status-queb
# or
PROFILE=queb ./scripts/dev status
```

## Profiles

The app supports domain-specific profiles:

- **Generic**: No domain configuration (localhost only)
- **Queb**: `config/profiles/queb.env` (queb.space domain)

### Using Profiles

1. **Environment Variable**: `PROFILE=queb ./scripts/dev serve`
2. **Command Line**: `./scripts/dev --profile=queb serve`
3. **Profile Scripts**: `./scripts/dev-queb serve`, `./scripts/deploy-queb`, `./scripts/status-queb`

## Configuration

- **Environment**: `.env` file for API keys
- **Profiles**: `config/profiles/[name].env` for domain settings
- **Project**: `config/project.js` for project metadata
- **SSL**: Auto-generated certificates for local development

## SSL Certificates

Generate trusted SSL certificates for local development:

```bash
# Generate certificates (stored in config/certs/)
./scripts/generate-certs.js

# Generate with custom domains
./scripts/generate-certs.js myapp.local dev.myapp.local

# Force regeneration
./scripts/generate-certs.js --force
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