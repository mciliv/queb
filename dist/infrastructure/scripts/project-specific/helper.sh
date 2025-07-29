#!/bin/sh

# Python Environment Setup for Molecular Analysis Project
# ======================================================
# This script uses the mol-toolkit utilities for Python environment management

set -e

# Load toolkit utilities
TOOLKIT_DIR="$(dirname "$(dirname "$(dirname "${BASH_SOURCE[0]}")")")/.mol-toolkit"
source "${TOOLKIT_DIR}/core/logging-utils.sh"
source "${TOOLKIT_DIR}/core/python-utils.sh"

# Load project configuration
source "$(dirname "$0")/../config.sh"

# Setup logging
setup_logging "$(basename "$0")"

log_info "Setting up Python environment for $PROJECT_NAME"
log_info "Python version: $PYTHON_VERSION"
log_info "Project directory: $PROJECT_ROOT"

# Setup Python environment using toolkit utilities
setup_python_env "$PYTHON_VERSION" "$PROJECT_ROOT"

log_info "Python environment setup completed successfully"