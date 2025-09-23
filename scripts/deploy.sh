#!/bin/bash

# Molecular Analysis App - GCP Deployment Script
# Deploys to Google Cloud Platform

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}🚀 Starting GCP deployment...${NC}"

# Check if we're in the right directory
if [ ! -f "package.json" ]; then
    echo -e "${RED}❌ Error: package.json not found. Are you in the project root?${NC}"
    exit 1
fi

# Check if gcloud is installed
if ! command -v gcloud &> /dev/null; then
    echo -e "${RED}❌ Error: gcloud CLI not found. Please install it first:${NC}"
    echo -e "${YELLOW}   https://cloud.google.com/sdk/docs/install${NC}"
    exit 1
fi

# Check if user is authenticated
if ! gcloud auth list --filter=status:ACTIVE --format="value(account)" | grep -q .; then
    echo -e "${RED}❌ Error: Not authenticated with GCP. Please run:${NC}"
    echo -e "${YELLOW}   gcloud auth login${NC}"
    exit 1
fi

# Get project ID
PROJECT_ID=$(gcloud config get-value project 2>/dev/null || echo "")
if [ -z "$PROJECT_ID" ]; then
    echo -e "${RED}❌ Error: No GCP project selected. Please run:${NC}"
    echo -e "${YELLOW}   gcloud config set project YOUR_PROJECT_ID${NC}"
    exit 1
fi

echo -e "${BLUE}📍 GCP Project: ${PROJECT_ID}${NC}"

# Build the project
echo -e "${BLUE}🔨 Building project...${NC}"
npm run build

# Check if build was successful
if [ ! -d "src/client/dist" ]; then
    echo -e "${RED}❌ Error: Build failed - dist directory not found${NC}"
    exit 1
fi

echo -e "${GREEN}✅ Build completed successfully${NC}"

# Deploy to GCP Cloud Functions
echo -e "${BLUE}☁️  Deploying to Google Cloud Functions...${NC}"

# Check if app.yaml exists for App Engine
if [ -f "app.yaml" ]; then
    echo -e "${BLUE}📱 Deploying to App Engine...${NC}"
    gcloud app deploy --quiet
    echo -e "${GREEN}✅ App Engine deployment completed${NC}"
    echo -e "${YELLOW}🌐 Your app is available at: https://${PROJECT_ID}.appspot.com${NC}"
else
    echo -e "${BLUE}⚡ Deploying to Cloud Functions...${NC}"
    
    # Deploy the main function
    gcloud functions deploy molecularAnalysis \
        --runtime nodejs18 \
        --trigger-http \
        --allow-unauthenticated \
        --source . \
        --entry-point molecularAnalysis \
        --memory 512MB \
        --timeout 60s \
        --quiet
    
    echo -e "${GREEN}✅ Cloud Functions deployment completed${NC}"
    
    # Get the function URL
    FUNCTION_URL=$(gcloud functions describe molecularAnalysis --format="value(httpsTrigger.url)" 2>/dev/null || echo "")
    if [ -n "$FUNCTION_URL" ]; then
        echo -e "${YELLOW}🌐 Your function is available at: ${FUNCTION_URL}${NC}"
    fi
fi

# Deploy static files to Cloud Storage (if needed)
echo -e "${BLUE}📦 Deploying static files to Cloud Storage...${NC}"

# Create a bucket for static files if it doesn't exist
BUCKET_NAME="${PROJECT_ID}-molecular-analysis-static"
if ! gsutil ls -b gs://${BUCKET_NAME} &> /dev/null; then
    echo -e "${BLUE}🪣 Creating Cloud Storage bucket: ${BUCKET_NAME}${NC}"
    gsutil mb gs://${BUCKET_NAME}
    gsutil iam ch allUsers:objectViewer gs://${BUCKET_NAME}
fi

# Upload static files
gsutil -m cp -r src/client/dist/* gs://${BUCKET_NAME}/
gsutil -m cp -r public/* gs://${BUCKET_NAME}/ 2>/dev/null || true

echo -e "${GREEN}✅ Static files uploaded to Cloud Storage${NC}"
echo -e "${YELLOW}🌐 Static files available at: https://storage.googleapis.com/${BUCKET_NAME}/index.html${NC}"

# Display deployment summary
echo -e "${GREEN}🎉 GCP deployment completed successfully!${NC}"
echo ""
echo -e "${BLUE}📋 Deployment Summary:${NC}"
echo -e "${YELLOW}   • Project: ${PROJECT_ID}${NC}"
if [ -f "app.yaml" ]; then
    echo -e "${YELLOW}   • App Engine: https://${PROJECT_ID}.appspot.com${NC}"
else
    echo -e "${YELLOW}   • Cloud Functions: ${FUNCTION_URL}${NC}"
fi
echo -e "${YELLOW}   • Static Files: https://storage.googleapis.com/${BUCKET_NAME}/index.html${NC}"
echo ""
echo -e "${BLUE}📄 To view logs: gcloud functions logs read molecularAnalysis${NC}"
echo -e "${BLUE}📄 To update: run 'ship' again${NC}"