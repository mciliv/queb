# GCloud Func - Google Cloud Functions Deployment Toolkit

A flat, specific toolkit for Google Cloud Functions deployment and management.

## Structure

```
gcloud-func/
├── deploy.sh          # Main deployment script
├── utils.sh           # All utility functions
├── config.sh          # Configuration template
├── auth.sh            # Authentication utilities
├── domain.sh          # Domain management
├── ssl.sh             # SSL certificate utilities
├── monitor.sh         # Monitoring and logging
└── README.md          # This file
```

## Quick Start

```bash
# Copy to your project
cp -r gcloud-func/ your-project/

# Update configuration
nano gcloud-func/config.sh

# Deploy
./gcloud-func/deploy.sh
```

## Usage

### Basic Deployment

```bash
#!/bin/bash

# Load utilities
source gcloud-func/utils.sh
source gcloud-func/config.sh

# Setup logging
setup_logging "$(basename "$0")"

# Deploy function
deploy_function "$FUNCTION_NAME" "$REGION" "$SOURCE_DIR"
```

### Domain Setup

```bash
#!/bin/bash

source gcloud-func/utils.sh
source gcloud-func/config.sh
source gcloud-func/domain.sh

# Setup domain
setup_domain "$DOMAIN_NAME"
check_dns_zone "$DNS_ZONE_NAME"
```

### SSL Certificate

```bash
#!/bin/bash

source gcloud-func/utils.sh
source gcloud-func/config.sh
source gcloud-func/ssl.sh

# Setup SSL
setup_ssl "$DOMAIN_NAME"
check_ssl "$DOMAIN_NAME"
```

## Available Functions

### Core Functions (utils.sh)
- `deploy_function name region source` - Deploy cloud function
- `check_function name region` - Check function status
- `get_function_url name region` - Get function URL
- `check_auth` - Check gcloud authentication
- `get_project` - Get current project ID
- `setup_logging script_name` - Setup logging
- `log_info message` - Log info message
- `log_error message` - Log error message

### Authentication (auth.sh)
- `login_gcloud` - Login to Google Cloud
- `set_project project_id` - Set active project
- `check_permissions` - Check required permissions
- `setup_service_account` - Setup service account

### Domain Management (domain.sh)
- `check_domain domain` - Check domain resolution
- `setup_domain domain` - Setup domain with DNS
- `check_dns_zone zone_name` - Check DNS zone
- `update_nameservers domain` - Update nameservers
- `test_url url` - Test URL accessibility

### SSL Certificates (ssl.sh)
- `setup_ssl domain` - Setup SSL certificate
- `check_ssl domain` - Check SSL status
- `renew_ssl domain` - Renew SSL certificate
- `force_ssl domain` - Force HTTPS redirect

### Monitoring (monitor.sh)
- `setup_monitoring function_name` - Setup function monitoring
- `check_logs function_name` - Check function logs
- `setup_alerts function_name` - Setup alerting
- `get_metrics function_name` - Get function metrics

## Configuration

Update `config.sh` with your settings:

```bash
# Project Configuration
PROJECT_NAME="your-project"
FUNCTION_NAME="your-function"
REGION="us-central1"
SOURCE_DIR="."

# Domain Configuration
DOMAIN_NAME="your-domain.com"
DNS_ZONE_NAME="your-dns-zone"

# Function Configuration
RUNTIME="nodejs20"
MEMORY="1GB"
TIMEOUT="540s"

# Export for use in scripts
export PROJECT_NAME FUNCTION_NAME REGION SOURCE_DIR DOMAIN_NAME DNS_ZONE_NAME RUNTIME MEMORY TIMEOUT
```

## Examples

### Deploy with Domain

```bash
#!/bin/bash

source gcloud-func/utils.sh
source gcloud-func/config.sh
source gcloud-func/domain.sh
source gcloud-func/ssl.sh

# Deploy function
deploy_function "$FUNCTION_NAME" "$REGION" "$SOURCE_DIR"

# Setup domain and SSL
setup_domain "$DOMAIN_NAME"
setup_ssl "$DOMAIN_NAME"

# Test deployment
test_url "https://$DOMAIN_NAME"
```

### Monitor Function

```bash
#!/bin/bash

source gcloud-func/utils.sh
source gcloud-func/config.sh
source gcloud-func/monitor.sh

# Setup monitoring
setup_monitoring "$FUNCTION_NAME"
setup_alerts "$FUNCTION_NAME"

# Check status
check_function "$FUNCTION_NAME" "$REGION"
check_logs "$FUNCTION_NAME"
```

## Benefits

- **Flat Structure**: Easy to copy and use
- **Specific**: Focused on Google Cloud Functions
- **Modular**: Load only what you need
- **Portable**: Works across projects
- **Well-documented**: Clear examples

## Contributing

When adding new utilities:
1. Add functions to appropriate .sh file
2. Include documentation in comments
3. Update this README with new functions
4. Add examples if applicable 