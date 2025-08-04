#!/bin/bash

# Domain Check Script for Molecular Analysis Project
# =================================================
# This script uses the mol-toolkit utilities for domain management

# Load toolkit utilities
TOOLKIT_DIR="$(dirname "$(dirname "$(dirname "${BASH_SOURCE[0]}")")")/.mol-toolkit"
source "${TOOLKIT_DIR}/core/logging-utils.sh"
source "${TOOLKIT_DIR}/core/domain-utils.sh"
source "${TOOLKIT_DIR}/core/gcloud-utils.sh"

# Load project configuration
source "$(dirname "$0")/../config.sh"

# Setup error handling and logging
set_error_handling "$(basename "$0")"

# Main script logic
log_info "Starting domain check for $DOMAIN_NAME"
log_info "Project: $PROJECT_NAME"
log_info "Region: $REGION"

echo "🌐 Domain Setup Status for $DOMAIN_NAME"
echo "======================================"
echo ""

# Check Google Cloud authentication
log_info "Checking Google Cloud authentication..."
if ! check_gcloud_auth; then
    log_error "Google Cloud not authenticated. Run: gcloud auth login"
    exit 1
fi

# Check DNS Zone Status
echo "📍 DNS Zone Status:"
check_gcloud_dns_zone "$DNS_ZONE_NAME" "$PROJECT_ID"
echo ""

# Check DNS Records
echo "📋 DNS Records:"
if gcloud dns record-sets list --zone="$DNS_ZONE_NAME" --project="$PROJECT_ID" --format="table(name,type,ttl,rrdatas)" 2>/dev/null; then
    echo "✅ DNS records configured"
else
    echo "❌ Cannot list DNS records"
fi
echo ""

# Check Domain Verification Status
echo "🔍 Domain Verification Status:"
check_domain_verification "$DOMAIN_NAME" "$PROJECT_ID"
echo ""

# Check Domain Mappings
echo "🔗 Domain Mappings:"
if ! check_domain_mappings "$FUNCTION_NAME" "$REGION" "$PROJECT_ID"; then
    echo "   → Run after domain verification:"
    echo "   → gcloud beta run domain-mappings create --service=$FUNCTION_NAME --domain=$DOMAIN_NAME --region=$REGION"
    echo "   → gcloud beta run domain-mappings create --service=$FUNCTION_NAME --domain=www.$DOMAIN_NAME --region=$REGION"
fi
echo ""

# Check Cloud Function Status
echo "☁️ Cloud Function Status:"
check_cloud_function "$FUNCTION_NAME" "$REGION" "$PROJECT_ID"
echo ""

# Check Nameserver Configuration
echo "🔧 Nameserver Configuration:"
echo "Current Google Cloud nameservers:"
get_gcloud_nameservers "$DNS_ZONE_NAME" "$PROJECT_ID" | sed 's/^/   /'
echo ""

# Check DNS Propagation
echo "🌍 DNS Propagation Check:"
check_domain_resolution "$DOMAIN_NAME"
check_domain_resolution "www.$DOMAIN_NAME"
echo ""

# Provide actionable next steps
echo "✅ Next Steps:"
echo "1. Update nameservers at Namecheap to Google Cloud nameservers above"
echo "2. Wait 24-48 hours for DNS propagation"
echo "3. Complete domain verification in Google Search Console"
echo "4. Run domain mapping commands shown above"
echo "5. Test: https://$DOMAIN_NAME and https://www.$DOMAIN_NAME"
echo ""

# Test current function URL
echo "🔗 Current Function URL Test:"
FUNCTION_URL=$(get_function_url "$FUNCTION_NAME" "$REGION" "$PROJECT_ID")
if [ -n "$FUNCTION_URL" ]; then
    test_url "$FUNCTION_URL"
else
    echo "⚠️ Cannot get function URL - checking function status..."
    gcloud functions describe "$FUNCTION_NAME" --region="$REGION" --project="$PROJECT_ID" --format="value(status)" 2>/dev/null || echo "❌ Function not found"
fi

log_info "Domain check completed" 