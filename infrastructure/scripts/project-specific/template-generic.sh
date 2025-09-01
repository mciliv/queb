#!/bin/bash
source "$(dirname "$0")/config.sh"
UTIL_DIR="/Users/m/util"
source "$UTIL_DIR/gcloud-func/utils.sh"
SCRIPT_NAME="$(basename "$0")"
LOG_FILE="${PROJECT_ROOT}/logs/${SCRIPT_NAME%.*}.log"
mkdir -p "$(dirname "$LOG_FILE")"
set -e
trap 'log_error "ERROR: Script failed at line $LINENO"' ERR
log_info "Starting deployment of $PROJECT_NAME to $DOMAIN_NAME"
log_info "Working directory: $PROJECT_ROOT"
log_info "Checking Google Cloud authentication..."
check_auth || {
    log_error "Please authenticate with Google Cloud: gcloud auth login"
    exit 1
}

log_info "Running pre-deployment checks..."
check_domain "$DOMAIN_NAME"
check_domain "www.$DOMAIN_NAME"

log_info "🚀 Starting deployment to Google Cloud Functions..."
log_info "Deploying with Cloud Functions Gen 1..."
# Load OPENAI_API_KEY from .env file if available
if [ -f "$PROJECT_ROOT/.env" ]; then
    log_info "📄 Loading environment variables from .env file..."
    OPENAI_API_KEY=$(grep -E "^OPENAI_API_KEY=" "$PROJECT_ROOT/.env" | cut -d'=' -f2- | tr -d '"' | tr -d "'")
    if [ -n "$OPENAI_API_KEY" ]; then
        log_info "✅ Found OPENAI_API_KEY in .env file"
    fi
fi

# Fallback to environment variable
OPENAI_API_KEY=${OPENAI_API_KEY:-""}
if [ -z "$OPENAI_API_KEY" ]; then
    log_warn "⚠️ OPENAI_API_KEY not found - molecular analysis features will not work"
    log_info "📝 API key should be in .env file or environment:"
    log_info "   1. Check .env file: OPENAI_API_KEY=your-key-here"
    log_info "   2. Or export locally: export OPENAI_API_KEY='your-key-here'"
    log_info "   3. Get key from: https://platform.openai.com/api-keys"
    log_info ""
    OPENAI_API_KEY="dummy-key-for-deployment"
else
    log_info "✅ OPENAI_API_KEY loaded successfully"
fi

gcloud functions deploy "$FUNCTION_NAME" \
    --runtime="nodejs20" \
    --source="$PROJECT_ROOT" \
    --entry-point="molecularAnalysis" \
    --trigger-http \
    --allow-unauthenticated \
    --memory="1GB" \
    --timeout="540s" \
    --region="$REGION" \
    --set-env-vars="FUNCTION_NAME=molecular-analysis,FUNCTION_TARGET=molecularAnalysis,NODE_ENV=production,OPENAI_API_KEY=$OPENAI_API_KEY"

log_info "Verifying deployment..."
url=$(gcloud functions describe "$FUNCTION_NAME" --region="$REGION" --format="value(url)" 2>/dev/null)
if [ -n "$url" ]; then
    log_info "✅ Function deployed successfully: $url"
    
    log_info "Testing deployed function..."
    test_url "$url"
    log_info "🌐 Verifying domain mapping for $DOMAIN_NAME..."
    run_url=$(gcloud run services describe "$FUNCTION_NAME" --region="$REGION" --format="value(status.url)" 2>/dev/null)
    
    if [ -n "$run_url" ]; then
        log_info "Testing Cloud Run service..."
        test_url "$run_url/health"
        
        log_info "Testing live domains..."
        if test_url "https://$DOMAIN_NAME/health"; then
            log_info "🎉 Deployment to $DOMAIN_NAME completed successfully!"
            log_info ""
            log_info "🌐 Your app is now live at:"
            log_info "  ✅ https://$DOMAIN_NAME"
            log_info "  ✅ https://www.$DOMAIN_NAME"
            log_info ""
            log_info "📡 Backend URLs:"
            log_info "  - $url (Cloud Function)"
            log_info "  - $run_url (Cloud Run)"
            log_info ""
            log_info "🔗 Try these endpoints:"
            log_info "  - https://$DOMAIN_NAME/ (Full App)"
            log_info "  - https://$DOMAIN_NAME/health (API status)"
            
            if [ "$OPENAI_API_KEY" = "dummy-key-for-deployment" ]; then
                log_info ""
                log_warn "🔑 Remember to set your OPENAI_API_KEY for full functionality:"
                log_info "   https://console.cloud.google.com/functions/details/us-central1/$FUNCTION_NAME?project=$PROJECT_ID&tab=variables"
            fi
        else
            log_warn "⚠️ Domain may need time to propagate. Check again in a few minutes."
            log_info "App deployed successfully at:"
            log_info "  - $url (direct Cloud Function)"
            log_info "  - $run_url (Cloud Run service)"
        fi
    else
        log_warn "⚠️ Could not get Cloud Run service URL, but deployment succeeded"
        log_info "Your app is live at: $url"
    fi
else
    log_error "❌ Function deployment verification failed"
    exit 1
fi

log_info "Deployment completed successfully" 