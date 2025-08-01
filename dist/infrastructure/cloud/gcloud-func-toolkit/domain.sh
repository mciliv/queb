#!/bin/bash

# GCloud Domain Management Utilities
# =================================

# Check domain resolution
check_domain() {
    local domain="$1"
    local timeout="${2:-5}"
    
    log_info "Checking domain: $domain"
    if command -v dig &> /dev/null; then
        if dig "$domain" A +short +timeout="$timeout" 2>/dev/null | grep -q .; then
            log_info "✅ Domain $domain is resolving"
            return 0
        else
            log_error "❌ Domain $domain is not resolving"
            return 1
        fi
    else
        log_warn "⚠️ dig command not available - install bind-utils"
        return 2
    fi
}

# Setup domain with DNS
setup_domain() {
    local domain="$1"
    
    log_info "Setting up domain: $domain"
    
    # Check if domain is verified
    local verified_domain=$(gcloud domains list-user-verified --project="$PROJECT_ID" --filter="domain:$domain" --format="value(domain)" 2>/dev/null)
    if [ -z "$verified_domain" ]; then
        log_info "Domain not verified. Please verify ownership in Google Cloud Console"
        log_info "Go to: https://console.cloud.google.com/apis/credentials/domainverification"
        return 1
    fi
    
    log_info "✅ Domain verified: $verified_domain"
    
    # Check DNS zone
    if ! check_dns_zone "$DNS_ZONE_NAME"; then
        log_info "Creating DNS zone: $DNS_ZONE_NAME"
        gcloud dns managed-zones create "$DNS_ZONE_NAME" \
            --dns-name="$domain." \
            --description="DNS zone for $domain" \
            --project="$PROJECT_ID"
    fi
    
    # Get nameservers
    local nameservers=$(get_gcloud_nameservers "$DNS_ZONE_NAME")
    if [ -n "$nameservers" ]; then
        log_info "Nameservers for $domain:"
        echo "$nameservers" | while read ns; do
            echo "   $ns"
        done
        log_warn "⚠️ Update your domain registrar with these nameservers"
    fi
    
    return 0
}

# Check DNS zone status
check_dns_zone() {
    local zone_name="$1"
    
    log_info "Checking DNS zone: $zone_name"
    if gcloud dns managed-zones describe "$zone_name" --project="$PROJECT_ID" --format="value(name,dnsName)" 2>/dev/null; then
        log_info "✅ DNS zone $zone_name exists"
        return 0
    else
        log_error "❌ DNS zone $zone_name not found"
        return 1
    fi
}

# Get Google Cloud nameservers
get_gcloud_nameservers() {
    local zone_name="$1"
    
    gcloud dns managed-zones describe "$zone_name" --project="$PROJECT_ID" --format="value(nameServers)" 2>/dev/null | tr ';' '\n'
}

# Test URL accessibility
test_url() {
    local url="$1"
    local timeout="${2:-5}"
    
    log_info "Testing URL: $url"
    if curl -s --max-time "$timeout" "$url" > /dev/null 2>&1; then
        log_info "✅ URL is responding"
        return 0
    else
        log_warn "❌ URL not responding (may be cold start)"
        return 1
    fi
} 