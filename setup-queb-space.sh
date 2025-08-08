#!/bin/bash

# Queb.Space Setup Script
# =======================
# Sets up both queb.space (Google Cloud) and dev.queb.space (local development)

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info() {
    echo -e "${BLUE}ℹ️  $1${NC}"
}

log_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

log_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

log_error() {
    echo -e "${RED}❌ $1${NC}"
}

# Get project root
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

echo "🚀 Queb.Space Setup"
echo "=================="
echo ""

# Check if .env exists
if [ ! -f ".env" ]; then
    log_warning "No .env file found. Creating template..."
    cat > .env << EOF
# Queb.Space Environment Configuration
NODE_ENV=development
PORT=8080

# OpenAI API (Required for molecular analysis)
OPENAI_API_KEY=your_openai_api_key_here

# Database Configuration
DB_HOST=localhost
DB_PORT=5432
DB_NAME=mol_users
DB_USER=mol_user
DB_PASSWORD=mol_password

# Google Cloud Configuration
GOOGLE_CLOUD_PROJECT=molecular-analysis
FUNCTION_NAME=molecular-analysis
FUNCTION_TARGET=molecularAnalysis

# Payment Configuration
PAYMENTS_ENABLED=false
PAYMENTS_DEV_MODE=true
PAYMENTS_REQUIRED=false

# Development Configuration
INTEGRATION_TEST=false
EOF
    log_success "Created .env template"
    log_warning "Please update OPENAI_API_KEY in .env file"
    echo ""
fi

# Function to setup local development (localhost)
setup_local_dev() {
    log_info "Setting up local development (localhost)..."
    
    # Check SSL certificates for localhost
    log_info "Checking SSL certificates for localhost..."
    if [ ! -d "backend/api/certs" ]; then
        mkdir -p backend/api/certs
    fi
    
    # Check if certificates already exist
    if [ ! -f "backend/api/certs/cert.pem" ] || [ ! -f "backend/api/certs/key.pem" ]; then
        log_info "SSL certificates will be generated automatically by the server"
        log_info "The server includes built-in HTTPS support"
    else
        log_success "SSL certificates already exist"
    fi
    
    # Hosts update no longer required; using localhost
    log_success "Hosts file changes not required"
    
    # Install missing dependencies if needed
    if ! command -v browser-sync >/dev/null 2>&1; then
        log_info "Installing browser-sync for live reload..."
        npm install -g browser-sync
    fi
    
    log_success "Local development setup complete"
    echo "   🌐 Access: http://localhost:3000 (or https://localhost:3002)"
    echo "   🔧 Start: ./dev-queb-space"
    echo ""
}

# Function to setup Google Cloud (queb.space)
setup_google_cloud() {
    log_info "Setting up Google Cloud (queb.space)..."
    
    # Check if gcloud is installed
    if ! command -v gcloud >/dev/null 2>&1; then
        log_error "Google Cloud SDK not found. Please install it first:"
        echo "   https://cloud.google.com/sdk/docs/install"
        return 1
    fi
    
    # Check authentication
    if ! gcloud auth list --filter=status:ACTIVE --format="value(account)" | grep -q .; then
        log_warning "Not authenticated with Google Cloud. Please run:"
        echo "   gcloud auth login"
        echo "   gcloud config set project molecular-analysis"
        return 1
    fi
    
    # Check if project is set
    CURRENT_PROJECT=$(gcloud config get-value project 2>/dev/null)
    if [ "$CURRENT_PROJECT" != "molecular-analysis" ]; then
        log_info "Setting project to molecular-analysis..."
        gcloud config set project molecular-analysis
    fi
    
    # Enable required APIs
    log_info "Enabling required Google Cloud APIs..."
    gcloud services enable cloudfunctions.googleapis.com
    gcloud services enable cloudbuild.googleapis.com
    gcloud services enable dns.googleapis.com
    gcloud services enable run.googleapis.com
    
    # Setup DNS zone for queb.space
    log_info "Setting up DNS zone for queb.space..."
    if ! gcloud dns managed-zones describe queb-space-zone --format="value(name)" 2>/dev/null; then
        gcloud dns managed-zones create queb-space-zone \
            --dns-name="queb.space." \
            --description="DNS zone for queb.space"
        log_success "Created DNS zone: queb-space-zone"
    else
        log_success "DNS zone already exists"
    fi
    
    # Get nameservers
    log_info "Google Cloud nameservers for queb.space:"
    gcloud dns managed-zones describe queb-space-zone --format="value(nameServers)" | tr ';' '\n' | sed 's/^/   /'
    echo ""
    log_warning "Update your domain registrar (Namecheap) with these nameservers"
    echo ""
    
    # Deploy the function
    log_info "Deploying to Google Cloud Functions..."
    if [ -z "$OPENAI_API_KEY" ]; then
        log_error "OPENAI_API_KEY not set. Please add it to .env file"
        return 1
    fi
    
    gcloud functions deploy molecular-analysis \
        --gen2 \
        --runtime=nodejs20 \
        --region=us-central1 \
        --source=. \
        --entry-point=molecularAnalysis \
        --trigger-http \
        --allow-unauthenticated \
        --memory=1GB \
        --timeout=540s \
        --set-env-vars="OPENAI_API_KEY=$OPENAI_API_KEY" \
        --quiet
    
    log_success "Function deployed successfully"
    
    # Get function URL
    FUNCTION_URL=$(gcloud functions describe molecular-analysis --region=us-central1 --format="value(serviceConfig.uri)")
    log_success "Function URL: $FUNCTION_URL"
    
    # Setup domain mapping
    log_info "Setting up domain mapping..."
    gcloud beta run domain-mappings create \
        --service=molecular-analysis \
        --domain=queb.space \
        --region=us-central1 \
        --quiet 2>/dev/null || log_warning "Domain mapping may already exist"
    
    gcloud beta run domain-mappings create \
        --service=molecular-analysis \
        --domain=www.queb.space \
        --region=us-central1 \
        --quiet 2>/dev/null || log_warning "Domain mapping may already exist"
    
    log_success "Google Cloud setup complete"
    echo "   🌐 Production: https://queb.space"
    echo "   🔧 Deploy: ./infrastructure/deployment/gcp"
    echo ""
}

# Main setup logic
echo "Choose setup option:"
echo "1. Local development only (dev.queb.space)"
echo "2. Google Cloud only (queb.space)"
echo "3. Both environments"
echo "4. Check status"
read -p "Enter choice (1-4): " choice

case $choice in
    1)
        setup_local_dev
        ;;
    2)
        setup_google_cloud
        ;;
    3)
        setup_local_dev
        echo ""
        setup_google_cloud
        ;;
    4)
        log_info "Checking setup status..."
        echo ""
        echo "🔍 Local Development Status:"
        if [ -f "backend/api/certs/cert.pem" ]; then
            log_success "SSL certificates: OK"
        else
            log_error "SSL certificates: Missing"
        fi
        
        if grep -q "dev.queb.space" /etc/hosts; then
            log_success "Hosts file: OK"
        else
            log_error "Hosts file: Missing entry"
        fi
        
        echo ""
        echo "☁️ Google Cloud Status:"
        if gcloud auth list --filter=status:ACTIVE --format="value(account)" | grep -q .; then
            log_success "Authentication: OK"
        else
            log_error "Authentication: Required"
        fi
        
        if gcloud dns managed-zones describe queb-space-zone --format="value(name)" 2>/dev/null; then
            log_success "DNS zone: OK"
        else
            log_error "DNS zone: Missing"
        fi
        
        if gcloud functions describe molecular-analysis --region=us-central1 --format="value(status)" 2>/dev/null; then
            log_success "Cloud Function: OK"
        else
            log_error "Cloud Function: Not deployed"
        fi
        ;;
    *)
        log_error "Invalid choice"
        exit 1
        ;;
esac

echo ""
log_success "Setup complete!"
echo ""
echo "📋 Next Steps:"
echo "1. Update OPENAI_API_KEY in .env file"
echo "2. For local dev: ./dev-simple"
echo "3. For production: ./infrastructure/deployment/gcp"
echo "4. Check domain status: ./infrastructure/scripts/project-specific/check-domain-status.sh" 