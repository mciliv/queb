#!/bin/bash

# Generic Script Template
# ======================
# This template shows how to create reusable scripts that work across projects

# Load configuration (required for all scripts)
source "$(dirname "$0")/config.sh"

# Script-specific configuration
# Override defaults if needed for this specific script
SCRIPT_NAME="$(basename "$0")"
LOG_FILE="${PROJECT_ROOT}/logs/${SCRIPT_NAME%.*}.log"

# Ensure log directory exists
mkdir -p "$(dirname "$LOG_FILE")"

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Error handling
set -e
trap 'log "ERROR: Script failed at line $LINENO"' ERR

# Main script logic
log "Starting $SCRIPT_NAME for project: $DOMAIN_NAME"
log "Working directory: $PROJECT_ROOT"

# Example: Generic domain check function
check_domain() {
    local domain="$1"
    log "Checking domain: $domain"
    
    if dig "$domain" A +short 2>/dev/null | grep -q .; then
        log "✅ Domain $domain is resolving"
        return 0
    else
        log "❌ Domain $domain is not resolving"
        return 1
    fi
}

# Example: Generic Google Cloud function check
check_function() {
    local function_name="$1"
    local region="$2"
    
    log "Checking function: $function_name in region: $region"
    
    if gcloud functions describe "$function_name" --region="$region" --format="value(status)" 2>/dev/null | grep -q "ACTIVE"; then
        log "✅ Function $function_name is active"
        return 0
    else
        log "❌ Function $function_name is not active"
        return 1
    fi
}

# Usage example
log "Running generic checks..."
check_domain "$DOMAIN_NAME"
check_domain "www.$DOMAIN_NAME"
check_function "$FUNCTION_NAME" "$REGION"

log "Script completed successfully" 