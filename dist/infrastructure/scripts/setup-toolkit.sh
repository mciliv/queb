#!/bin/bash

# Setup Toolkit Script
# ===================
# Automatically sets up the mol-toolkit as a git submodule

# Load toolkit utilities (if they exist)
TOOLKIT_DIR="$(dirname "$0")/../.mol-toolkit"
if [ -d "$TOOLKIT_DIR" ]; then
    source "${TOOLKIT_DIR}/core/logging-utils.sh"
    set_error_handling "$(basename "$0")"
else
    # Fallback logging if toolkit doesn't exist yet
    log() {
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    }
    log_info() { log "INFO: $1"; }
    log_error() { log "ERROR: $1"; }
    log_warn() { log "WARN: $1"; }
fi

# Configuration
TOOLKIT_REPO="https://github.com/mciliv/mol-toolkit.git"
TOOLKIT_BRANCH="main"
TOOLKIT_DIR_NAME=".mol-toolkit"
PROJECT_ROOT="$(dirname "$(dirname "$0")")"

log_info "Setting up mol-toolkit for project: $(basename "$PROJECT_ROOT")"

# Check if we're in a git repository
if ! git rev-parse --git-dir > /dev/null 2>&1; then
    log_error "Not in a git repository. Please run this from your project root."
    exit 1
fi

# Check if toolkit already exists
if [ -d "$TOOLKIT_DIR_NAME" ]; then
    log_warn "Toolkit directory already exists: $TOOLKIT_DIR_NAME"
    echo "Options:"
    echo "1. Remove existing and reinstall"
    echo "2. Update existing toolkit"
    echo "3. Skip setup"
    read -p "Choose option (1-3): " choice
    
    case $choice in
        1)
            log_info "Removing existing toolkit..."
            rm -rf "$TOOLKIT_DIR_NAME"
            git rm -rf "$TOOLKIT_DIR_NAME" 2>/dev/null || true
            ;;
        2)
            log_info "Updating existing toolkit..."
            cd "$TOOLKIT_DIR_NAME"
            git pull origin "$TOOLKIT_BRANCH"
            cd "$PROJECT_ROOT"
            log_info "Toolkit updated successfully"
            exit 0
            ;;
        3)
            log_info "Skipping toolkit setup"
            exit 0
            ;;
        *)
            log_error "Invalid choice"
            exit 1
            ;;
    esac
fi

# Add toolkit as submodule
log_info "Adding mol-toolkit as git submodule..."
if git submodule add -b "$TOOLKIT_BRANCH" "$TOOLKIT_REPO" "$TOOLKIT_DIR_NAME"; then
    log_info "Submodule added successfully"
else
    log_error "Failed to add submodule"
    exit 1
fi

# Initialize and update submodule
log_info "Initializing submodule..."
git submodule update --init --recursive

# Make toolkit scripts executable
log_info "Setting up toolkit permissions..."
find "$TOOLKIT_DIR_NAME" -name "*.sh" -type f -exec chmod +x {} \;

# Create project-specific config if it doesn't exist
PROJECT_CONFIG="$(dirname "$0")/config.sh"
if [ ! -f "$PROJECT_CONFIG" ]; then
    log_info "Creating project configuration..."
    cat > "$PROJECT_CONFIG" << 'EOF'
#!/bin/bash

# Project Configuration
# ====================

# Project Information
PROJECT_NAME="$(basename "$(dirname "$(dirname "$0")")")"
PROJECT_DESCRIPTION="Your project description"

# Domain Configuration
DOMAIN_NAME="your-domain.com"
DNS_ZONE_NAME="your-dns-zone"
REGION="us-central1"
FUNCTION_NAME="your-function-name"

# Google Cloud Configuration
PROJECT_ID="${GOOGLE_CLOUD_PROJECT:-$(gcloud config get-value project 2>/dev/null)}"

# Python Configuration
PYTHON_VERSION="3.12.9"

# Paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Export for use in other scripts
export PROJECT_NAME PROJECT_DESCRIPTION DOMAIN_NAME DNS_ZONE_NAME REGION FUNCTION_NAME PROJECT_ID PYTHON_VERSION SCRIPT_DIR PROJECT_ROOT
EOF
    chmod +x "$PROJECT_CONFIG"
    log_info "Created project config: $PROJECT_CONFIG"
    log_warn "Please update the configuration with your project-specific values"
fi

# Create project-specific scripts directory if it doesn't exist
PROJECT_SCRIPTS_DIR="$(dirname "$0")/project-specific"
if [ ! -d "$PROJECT_SCRIPTS_DIR" ]; then
    log_info "Creating project-specific scripts directory..."
    mkdir -p "$PROJECT_SCRIPTS_DIR"
fi

# Test the setup
log_info "Testing toolkit setup..."
if [ -f "$TOOLKIT_DIR_NAME/core/logging-utils.sh" ]; then
    log_info "✅ Toolkit setup successful!"
    echo ""
    echo "Next steps:"
    echo "1. Update $PROJECT_CONFIG with your project-specific values"
    echo "2. Create project-specific scripts in $PROJECT_SCRIPTS_DIR"
    echo "3. Use toolkit functions in your scripts:"
    echo "   source \"$TOOLKIT_DIR_NAME/core/logging-utils.sh\""
    echo ""
    echo "Example:"
    echo "cp $TOOLKIT_DIR_NAME/templates/domain-check.sh.template $PROJECT_SCRIPTS_DIR/check-domain.sh"
else
    log_error "❌ Toolkit setup failed - logging utilities not found"
    exit 1
fi

# Commit the changes
log_info "Committing toolkit setup..."
git add "$TOOLKIT_DIR_NAME" "$PROJECT_CONFIG" "$PROJECT_SCRIPTS_DIR"
git commit -m "Add mol-toolkit as submodule and setup project configuration"

log_info "Toolkit setup completed successfully!" 