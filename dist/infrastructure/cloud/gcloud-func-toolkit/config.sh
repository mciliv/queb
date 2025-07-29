#!/bin/bash

# GCloud Functions Configuration for Mol Project
# =============================================

# Project Configuration
PROJECT_NAME="molecular-analysis"
FUNCTION_NAME="molecular-analysis"
REGION="us-central1"
SOURCE_DIR="."

# Domain Configuration
DOMAIN_NAME="queb.space"
DNS_ZONE_NAME="queb-space-zone"

# Function Configuration
RUNTIME="nodejs20"
MEMORY="1GB"
TIMEOUT="540s"
ENTRY_POINT="main"

# Google Cloud Configuration
PROJECT_ID="${GOOGLE_CLOUD_PROJECT:-$(gcloud config get-value project 2>/dev/null)}"

# Paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Export for use in scripts
export PROJECT_NAME FUNCTION_NAME REGION SOURCE_DIR DOMAIN_NAME DNS_ZONE_NAME RUNTIME MEMORY TIMEOUT ENTRY_POINT PROJECT_ID SCRIPT_DIR PROJECT_ROOT 