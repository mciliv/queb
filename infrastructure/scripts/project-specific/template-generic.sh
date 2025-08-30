#!/bin/bash

# Molecular Analysis Deployment Script
# ====================================
# Deploys the molecular analysis application to Google Cloud Functions

# Load configuration (required for all scripts)
source "$(dirname "$0")/config.sh"

# Load GCloud utilities from util directory
UTIL_DIR="/Users/m/util"
source "$UTIL_DIR/gcloud-func/utils.sh"

# Script-specific configuration
SCRIPT_NAME="$(basename "$0")"
LOG_FILE="${PROJECT_ROOT}/logs/${SCRIPT_NAME%.*}.log"

# Ensure log directory exists
mkdir -p "$(dirname "$LOG_FILE")"

# Error handling
set -e
trap 'log_error "ERROR: Script failed at line $LINENO"' ERR

# Main deployment logic
log_info "Starting deployment of $PROJECT_NAME to $DOMAIN_NAME"
log_info "Working directory: $PROJECT_ROOT"

# Check authentication
log_info "Checking Google Cloud authentication..."
check_auth || {
    log_error "Please authenticate with Google Cloud: gcloud auth login"
    exit 1
}

# Check domain resolution
log_info "Running pre-deployment checks..."
check_domain "$DOMAIN_NAME"
check_domain "www.$DOMAIN_NAME"

# Deploy the function
log_info "🚀 Starting deployment to Google Cloud Functions..."

# Deploy with functions.js entry point and environment variables
gcloud functions deploy "$FUNCTION_NAME" \
    --region="$REGION" \
    --source="$PROJECT_ROOT" \
    --runtime="nodejs20" \
    --memory="1GB" \
    --timeout="540s" \
    --trigger-http \
    --allow-unauthenticated \
    --entry-point="molecularAnalysis" \
    --set-env-vars="FUNCTION_NAME=molecular-analysis,FUNCTION_TARGET=molecularAnalysis,NODE_ENV=production"

# Verify deployment
log_info "Verifying deployment..."
if check_function "$FUNCTION_NAME" "$REGION"; then
    local url=$(gcloud functions describe "$FUNCTION_NAME" --region="$REGION" --format="value(httpsTrigger.url)" 2>/dev/null)
    log_info "✅ Function deployed successfully: $url"
    
    # Test the deployed function
    log_info "Testing deployed function..."
    test_url "$url"
    
    log_info "🎉 Deployment to $DOMAIN_NAME completed successfully!"
else
    log_error "❌ Function deployment verification failed"
    exit 1
fi

log_info "Deployment completed successfully" 