#!/bin/bash

# Simple Project Configuration
# ============================

# Project paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Export for use in scripts
export SCRIPT_DIR PROJECT_ROOT 