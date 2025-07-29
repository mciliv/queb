#!/bin/bash

# GCloud SSL Certificate Utilities
# ===============================

# Setup SSL certificate
setup_ssl() {
    local domain="$1"
    
    log_info "Setting up SSL certificate for: $domain"
    
    # Check if certificate already exists
    local existing_cert=$(gcloud compute ssl-certificates list --project="$PROJECT_ID" --filter="name:$domain" --format="value(name)" 2>/dev/null)
    
    if [ -n "$existing_cert" ]; then
        log_info "✅ SSL certificate already exists: $existing_cert"
        return 0
    fi
    
    # Create SSL certificate
    log_info "Creating SSL certificate..."
    if gcloud compute ssl-certificates create "$domain" \
        --domains="$domain" \
        --project="$PROJECT_ID"; then
        log_info "✅ SSL certificate created: $domain"
        return 0
    else
        log_error "❌ Failed to create SSL certificate"
        return 1
    fi
}

# Check SSL certificate status
check_ssl() {
    local domain="$1"
    
    log_info "Checking SSL certificate for: $domain"
    
    local cert_status=$(gcloud compute ssl-certificates describe "$domain" --project="$PROJECT_ID" --format="value(managed.status)" 2>/dev/null)
    
    if [ -n "$cert_status" ]; then
        case "$cert_status" in
            "ACTIVE")
                log_info "✅ SSL certificate is active"
                return 0
                ;;
            "PROVISIONING")
                log_info "⏳ SSL certificate is provisioning"
                return 0
                ;;
            "FAILED_NOT_VISIBLE")
                log_error "❌ SSL certificate failed - domain not visible"
                return 1
                ;;
            "FAILED")
                log_error "❌ SSL certificate failed"
                return 1
                ;;
            *)
                log_warn "⚠️ SSL certificate status: $cert_status"
                return 1
                ;;
        esac
    else
        log_error "❌ SSL certificate not found"
        return 1
    fi
}

# Wait for SSL certificate to be active
wait_for_ssl() {
    local domain="$1"
    local max_attempts="${2:-30}"
    local delay="${3:-10}"
    
    log_info "Waiting for SSL certificate to be active: $domain"
    
    for ((i=1; i<=max_attempts; i++)); do
        log_info "Attempt $i/$max_attempts"
        
        if check_ssl "$domain"; then
            log_info "✅ SSL certificate is active"
            return 0
        fi
        
        if [ $i -lt $max_attempts ]; then
            log_info "Waiting $delay seconds before next check..."
            sleep "$delay"
        fi
    done
    
    log_error "❌ SSL certificate activation timeout"
    return 1
} 