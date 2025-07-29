#!/bin/bash

# Update Toolkit Script
# ====================
# Updates the mol-toolkit submodule to the latest version

# Load toolkit utilities
TOOLKIT_DIR="$(dirname "$0")/../.mol-toolkit"
if [ -d "$TOOLKIT_DIR" ]; then
    source "${TOOLKIT_DIR}/core/logging-utils.sh"
    set_error_handling "$(basename "$0")"
else
    # Fallback logging if toolkit doesn't exist
    log() {
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    }
    log_info() { log "INFO: $1"; }
    log_error() { log "ERROR: $1"; }
    log_warn() { log "WARN: $1"; }
fi

TOOLKIT_DIR_NAME=".mol-toolkit"
PROJECT_ROOT="$(dirname "$(dirname "$0")")"

log_info "Updating mol-toolkit..."

# Check if we're in a git repository
if ! git rev-parse --git-dir > /dev/null 2>&1; then
    log_error "Not in a git repository. Please run this from your project root."
    exit 1
fi

# Check if toolkit exists
if [ ! -d "$TOOLKIT_DIR_NAME" ]; then
    log_error "Toolkit not found. Run setup-toolkit.sh first."
    exit 1
fi

# Get current commit hash
cd "$TOOLKIT_DIR_NAME"
CURRENT_COMMIT=$(git rev-parse HEAD)
CURRENT_BRANCH=$(git branch --show-current)
log_info "Current toolkit commit: $CURRENT_COMMIT ($CURRENT_BRANCH)"

# Fetch latest changes
log_info "Fetching latest changes..."
git fetch origin

# Check if there are updates
LATEST_COMMIT=$(git rev-parse origin/$CURRENT_BRANCH)
if [ "$CURRENT_COMMIT" = "$LATEST_COMMIT" ]; then
    log_info "✅ Toolkit is already up to date"
    exit 0
fi

# Show what's changing
log_info "Updates available:"
git log --oneline "$CURRENT_COMMIT..origin/$CURRENT_BRANCH"

# Ask for confirmation
echo ""
read -p "Update toolkit to latest version? (y/N): " confirm
if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
    log_info "Update cancelled"
    exit 0
fi

# Update the submodule
log_info "Updating toolkit..."
git pull origin "$CURRENT_BRANCH"

# Go back to project root
cd "$PROJECT_ROOT"

# Update submodule reference
log_info "Updating submodule reference..."
git add "$TOOLKIT_DIR_NAME"
git commit -m "Update mol-toolkit to latest version"

# Make scripts executable
log_info "Setting up permissions..."
find "$TOOLKIT_DIR_NAME" -name "*.sh" -type f -exec chmod +x {} \;

log_info "✅ Toolkit updated successfully!"
log_info "New commit: $(cd "$TOOLKIT_DIR_NAME" && git rev-parse HEAD)" 