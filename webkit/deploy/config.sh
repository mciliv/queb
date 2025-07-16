#!/bin/bash

# Simple App Deployment Toolkit Configuration
# ==========================================

# Project Information
PROJECT_NAME="your-project"
PROJECT_DESCRIPTION="Your project description"

# Domain Configuration
DOMAIN_NAME="your-domain.com"
DNS_ZONE_NAME="your-dns-zone"
REGION="us-central1"
FUNCTION_NAME="your-function-name"

# Google Cloud Configuration
PROJECT_ID="${GOOGLE_CLOUD_PROJECT:-$(gcloud config get-value project 2>/dev/null)}"

# Python Configuration
PYTHON_VERSION="3.12"

# Paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Export for use in other scripts
export PROJECT_NAME PROJECT_DESCRIPTION DOMAIN_NAME DNS_ZONE_NAME REGION FUNCTION_NAME PROJECT_ID PYTHON_VERSION SCRIPT_DIR PROJECT_ROOT 