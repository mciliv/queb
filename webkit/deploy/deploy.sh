#!/bin/bash

# Simple App Deployment Script
# ===========================

# Load toolkit
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/utils.sh"
source "$SCRIPT_DIR/config.sh"

# Setup logging
setup_logging "$(basename "$0")"

log_info "Starting deployment for $PROJECT_NAME"

# Validate environment
validate_env "PROJECT_NAME" "DOMAIN_NAME" "REGION" "FUNCTION_NAME" || {
    log_error "Missing required configuration. Please update config.sh"
    exit 1
}

# Check authentication
check_auth || {
    log_error "Please authenticate with Google Cloud: gcloud auth login"
    exit 1
}

# Check domain
check_domain "$DOMAIN_NAME"

# Deploy function
if [ -d "$PROJECT_ROOT" ]; then
    deploy_function "$FUNCTION_NAME" "$REGION" "$PROJECT_ROOT"
else
    log_error "Project root not found: $PROJECT_ROOT"
    exit 1
fi

# Test deployment
if [ $? -eq 0 ]; then
    local url=$(gcloud functions describe "$FUNCTION_NAME" --region="$REGION" --project="$PROJECT_ID" --format="value(httpsTrigger.url)" 2>/dev/null)
    if [ -n "$url" ]; then
        log_info "Testing deployed function..."
        test_url "$url"
    fi
fi

log_info "Deployment completed!" 