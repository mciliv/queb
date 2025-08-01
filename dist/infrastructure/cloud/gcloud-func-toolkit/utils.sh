#!/bin/bash

# Simple App Deployment Toolkit - All utilities in one file
# ========================================================

# Logging Functions
# ================

LOG_LEVEL=${LOG_LEVEL:-1}  # 0=debug, 1=info, 2=warn, 3=error

log() {
    local level="$1"
    local message="$2"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    local script_name=$(basename "${BASH_SOURCE[1]:-unknown}")
    
    case "$level" in
        "DEBUG") echo -e "\033[36m[${timestamp}] [DEBUG] [${script_name}] ${message}\033[0m" ;;
        "INFO")  echo -e "\033[32m[${timestamp}] [INFO]  [${script_name}] ${message}\033[0m" ;;
        "WARN")  echo -e "\033[33m[${timestamp}] [WARN]  [${script_name}] ${message}\033[0m" ;;
        "ERROR") echo -e "\033[31m[${timestamp}] [ERROR] [${script_name}] ${message}\033[0m" ;;
    esac
}

log_debug() { [ $LOG_LEVEL -le 0 ] && log "DEBUG" "$1"; }
log_info()  { [ $LOG_LEVEL -le 1 ] && log "INFO"  "$1"; }
log_warn()  { [ $LOG_LEVEL -le 2 ] && log "WARN"  "$1"; }
log_error() { [ $LOG_LEVEL -le 3 ] && log "ERROR" "$1"; }

setup_logging() {
    local script_name="$1"
    log_info "Logging initialized for $script_name"
}

# Domain & DNS Functions
# =====================

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

check_dns_zone() {
    local zone_name="$1"
    local project_id="${2:-$(gcloud config get-value project 2>/dev/null)}"
    
    log_info "Checking DNS zone: $zone_name"
    if gcloud dns managed-zones describe "$zone_name" --project="$project_id" --format="value(name,dnsName)" 2>/dev/null; then
        log_info "✅ DNS zone $zone_name exists"
        return 0
    else
        log_error "❌ DNS zone $zone_name not found"
        return 1
    fi
}

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

# Google Cloud Functions
# =====================

check_auth() {
    log_info "Checking Google Cloud authentication"
    if gcloud auth list --filter=status:ACTIVE --format="value(account)" 2>/dev/null | grep -q .; then
        log_info "✅ Google Cloud authenticated"
        return 0
    else
        log_error "❌ Google Cloud not authenticated"
        return 1
    fi
}

get_project() {
    gcloud config get-value project 2>/dev/null
}

check_function() {
    local function_name="$1"
    local region="$2"
    local project_id="${3:-$(get_project)}"
    
    log_info "Checking function: $function_name in $region"
    local function_status=$(gcloud functions describe "$function_name" --region="$region" --project="$project_id" --format="value(name,status)" 2>/dev/null)
    if [ -n "$function_status" ]; then
        log_info "✅ Cloud Function: $function_status"
        return 0
    else
        log_error "❌ Cloud Function $function_name not found"
        return 1
    fi
}

get_function_url() {
    local function_name="$1"
    local region="$2"
    local project_id="${3:-$(get_project)}"
    
    log_info "Getting function URL: $function_name in $region"
    local url=$(gcloud functions describe "$function_name" --region="$region" --project="$project_id" --format="value(httpsTrigger.url)" 2>/dev/null)
    if [ -n "$url" ]; then
        log_info "✅ Function URL: $url"
        echo "$url"
        return 0
    else
        log_error "❌ Could not get function URL"
        return 1
    fi
}

deploy_function() {
    local function_name="$1"
    local region="$2"
    local source_dir="$3"
    local project_id="${4:-$(get_project)}"
    local runtime="${5:-nodejs20}"
    local memory="${6:-1GB}"
    local timeout="${7:-540s}"
    local entry_point="${8:-main}"
    
    log_info "🚀 Deploying $function_name to $region with entry point: $entry_point..."
    
    gcloud functions deploy "$function_name" \
        --region="$region" \
        --source="$source_dir" \
        --runtime="$runtime" \
        --memory="$memory" \
        --timeout="$timeout" \
        --entry-point="$entry_point" \
        --trigger-http \
        --allow-unauthenticated \
        --project="$project_id"
    
    if [ $? -eq 0 ]; then
        local url=$(gcloud functions describe "$function_name" --region="$region" --project="$project_id" --format="value(httpsTrigger.url)" 2>/dev/null)
        log_info "✅ Deployed successfully: $url"
        return 0
    else
        log_error "❌ Deployment failed"
        return 1
    fi
}

# Python Environment Functions
# ===========================

check_python() {
    local version="$1"
    log_info "Checking Python version: $version"
    
    if command -v python3 &> /dev/null; then
        local current_version=$(python3 --version 2>&1 | cut -d' ' -f2)
        log_info "✅ Python $current_version is available"
        return 0
    else
        log_error "❌ Python not found"
        return 1
    fi
}

setup_python_env() {
    local version="$1"
    log_info "Setting up Python environment for version: $version"
    
    # Check if pyenv is available
    if command -v pyenv &> /dev/null; then
        log_info "Using pyenv to manage Python version"
        pyenv install -s "$version" 2>/dev/null || true
        pyenv local "$version"
    else
        log_warn "pyenv not available, using system Python"
    fi
    
    # Check if poetry is available
    if command -v poetry &> /dev/null; then
        log_info "✅ Poetry is available for dependency management"
    else
        log_warn "⚠️ Poetry not found - install with: curl -sSL https://install.python-poetry.org | python3 -"
    fi
    
    log_info "✅ Python environment setup complete"
}

install_deps() {
    local requirements_file="$1"
    
    log_info "Installing dependencies from: $requirements_file"
    
    if [ -f "$requirements_file" ]; then
        if command -v poetry &> /dev/null; then
            poetry install
        else
            pip install -r "$requirements_file"
        fi
        log_info "✅ Dependencies installed"
        return 0
    else
        log_error "❌ Requirements file not found: $requirements_file"
        return 1
    fi
}

# Utility Functions
# ================

validate_env() {
    local required_vars=("$@")
    local missing_vars=()
    
    for var in "${required_vars[@]}"; do
        if [ -z "${!var}" ]; then
            missing_vars+=("$var")
        fi
    done
    
    if [ ${#missing_vars[@]} -gt 0 ]; then
        log_error "❌ Missing required environment variables: ${missing_vars[*]}"
        return 1
    else
        log_info "✅ All required environment variables are set"
        return 0
    fi
}

show_progress() {
    local message="$1"
    local current="$2"
    local total="$3"
    
    local percentage=$((current * 100 / total))
    local filled=$((percentage / 2))
    local empty=$((50 - filled))
    
    local bar=""
    for ((i=0; i<filled; i++)); do bar+="█"; done
    for ((i=0; i<empty; i++)); do bar+="░"; done
    
    printf "\r%s [%s] %d%% (%d/%d)" "$message" "$bar" "$percentage" "$current" "$total"
    
    if [ "$current" -eq "$total" ]; then
        echo ""
    fi
} 