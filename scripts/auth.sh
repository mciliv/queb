#!/bin/bash

# GCloud Authentication Utilities
# ==============================

# Login to Google Cloud
login_gcloud() {
    log_info "Logging into Google Cloud..."
    if gcloud auth login --no-launch-browser; then
        log_info "✅ Successfully logged into Google Cloud"
        return 0
    else
        log_error "❌ Failed to login to Google Cloud"
        return 1
    fi
}

# Set active project
set_project() {
    local project_id="$1"
    log_info "Setting active project to: $project_id"
    if gcloud config set project "$project_id"; then
        log_info "✅ Project set to: $project_id"
        return 0
    else
        log_error "❌ Failed to set project: $project_id"
        return 1
    fi
}

# Check required permissions
check_permissions() {
    log_info "Checking required permissions..."
    
    local required_roles=(
        "roles/cloudfunctions.developer"
        "roles/iam.serviceAccountUser"
        "roles/dns.admin"
        "roles/cloudbuild.builds.builder"
    )
    
    local missing_roles=()
    
    for role in "${required_roles[@]}"; do
        if ! gcloud projects get-iam-policy "$PROJECT_ID" --flatten="bindings[].members" --format="table(bindings.role)" --filter="bindings.role:$role" | grep -q "$role"; then
            missing_roles+=("$role")
        fi
    done
    
    if [ ${#missing_roles[@]} -gt 0 ]; then
        log_warn "⚠️ Missing required roles:"
        for role in "${missing_roles[@]}"; do
            echo "   - $role"
        done
        log_info "Run: gcloud projects add-iam-policy-binding $PROJECT_ID --member=user:$(gcloud config get-value account) --role=$role"
        return 1
    else
        log_info "✅ All required permissions are set"
        return 0
    fi
}

# Setup service account
setup_service_account() {
    local service_account_name="$1"
    local service_account_email="$service_account_name@$PROJECT_ID.iam.gserviceaccount.com"
    
    log_info "Setting up service account: $service_account_name"
    
    # Create service account if it doesn't exist
    if ! gcloud iam service-accounts describe "$service_account_email" --project="$PROJECT_ID" >/dev/null 2>&1; then
        log_info "Creating service account: $service_account_name"
        gcloud iam service-accounts create "$service_account_name" \
            --display-name="$service_account_name" \
            --project="$PROJECT_ID"
    fi
    
    # Grant required roles
    local roles=(
        "roles/cloudfunctions.developer"
        "roles/iam.serviceAccountUser"
        "roles/dns.admin"
    )
    
    for role in "${roles[@]}"; do
        log_info "Granting role: $role"
        gcloud projects add-iam-policy-binding "$PROJECT_ID" \
            --member="serviceAccount:$service_account_email" \
            --role="$role"
    done
    
    # Create and download key
    local key_file="key-$service_account-name.json"
    log_info "Creating service account key: $key_file"
    gcloud iam service-accounts keys create "$key_file" \
        --iam-account="$service_account_email" \
        --project="$PROJECT_ID"
    
    log_info "✅ Service account setup complete"
    log_warn "⚠️ Keep the key file secure: $key_file"
}

# Verify authentication and project
verify_setup() {
    log_info "Verifying Google Cloud setup..."
    
    # Check authentication
    if ! check_auth; then
        log_error "❌ Not authenticated. Run: login_gcloud"
        return 1
    fi
    
    # Check project
    local current_project=$(get_project)
    if [ -z "$current_project" ]; then
        log_error "❌ No project set. Run: set_project PROJECT_ID"
        return 1
    fi
    
    log_info "✅ Using project: $current_project"
    
    # Check permissions
    if ! check_permissions; then
        log_warn "⚠️ Some permissions may be missing"
    fi
    
    return 0
} 