# Simple App Deployment Toolkit

A minimal, flat toolkit for basic app deployment operations.

## Structure

```
toolkit/
├── utils.sh          # All utility functions in one file
├── config.sh         # Configuration template
├── deploy.sh         # Basic deployment script
└── README.md         # This file
```

## Usage

### Quick Start

1. Copy the `toolkit/` directory to your project
2. Update `config.sh` with your project settings
3. Use `deploy.sh` or source `utils.sh` in your scripts

### In Your Scripts

```bash
#!/bin/bash

# Load toolkit
source "$(dirname "$0")/toolkit/utils.sh"
source "$(dirname "$0")/toolkit/config.sh"

# Setup logging
setup_logging "$(basename "$0")"

# Use functions
log_info "Starting deployment..."
check_domain "$DOMAIN_NAME"
deploy_function "$FUNCTION_NAME" "$REGION"
```

## Available Functions

### Domain & DNS
- `check_domain domain` - Check if domain resolves
- `check_dns_zone zone_name` - Check DNS zone status
- `test_url url` - Test URL accessibility

### Google Cloud
- `check_auth` - Check gcloud authentication
- `get_project` - Get current project ID
- `check_function name region` - Check function status
- `deploy_function name region source` - Deploy cloud function

### Python Environment
- `check_python version` - Check Python version
- `setup_python_env version` - Setup Python environment
- `install_deps requirements_file` - Install dependencies

### Logging
- `log_info message` - Log info message
- `log_error message` - Log error message
- `setup_logging script_name` - Setup logging

## Configuration

Update `config.sh` with your project settings:

```bash
PROJECT_NAME="your-project"
DOMAIN_NAME="your-domain.com"
REGION="us-central1"
FUNCTION_NAME="your-function"
PYTHON_VERSION="3.12"
```

## Benefits

- **Simple**: Single file for all utilities
- **Portable**: Easy to copy between projects
- **Minimal**: Only essential functions
- **Self-contained**: No complex dependencies 