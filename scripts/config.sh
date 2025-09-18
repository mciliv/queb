#!/bin/bash

# GCloud Functions Configuration
# =============================

# Project Configuration
PROJECT_NAME="your-project"
FUNCTION_NAME="your-function"
REGION="us-central1"
SOURCE_DIR="."

# Domain Configuration
DOMAIN_NAME="your-domain.com"
DNS_ZONE_NAME="your-dns-zone"

# Function Configuration
RUNTIME="nodejs20"
MEMORY="2GB"
TIMEOUT="540s"
ENTRY_POINT="molecularAnalysis"
# Use Dockerfile for more reliable deployment
SOURCE_DIR="."
DOCKERFILE="Dockerfile"

# Google Cloud Configuration
PROJECT_ID="${GOOGLE_CLOUD_PROJECT:-$(gcloud config get-value project 2>/dev/null)}"

# Paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Export for use in scripts
export PROJECT_NAME FUNCTION_NAME REGION SOURCE_DIR DOMAIN_NAME DNS_ZONE_NAME RUNTIME MEMORY TIMEOUT ENTRY_POINT PROJECT_ID SCRIPT_DIR PROJECT_ROOT 