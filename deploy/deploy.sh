#!/bin/bash

# Simple GCP Deployment Script
# Deploy to Google Cloud Platform
# Usage: ./deploy.sh [dev|prod]

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Determine environment
ENV=${1:-prod}
if [ "$ENV" = "dev" ]; then
  echo -e "${BLUE}🚀 Deploying to DEV environment...${NC}"
else
  echo -e "${BLUE}🚀 Deploying to PRODUCTION...${NC}"
fi

# Fast deployment - skip unnecessary checks
PROJECT_ID="mol-analysis-app"

# Build in parallel with deployment prep
echo -e "${BLUE}🔨 Building...${NC}"
if [ "$ENV" = "dev" ]; then
  # Development build with source maps and unminified code
  NODE_ENV=development node src/client/build-frontend.js --dev &
else
  # Production build
  node src/client/build-frontend.js &
fi
BUILD_PID=$!

# Prepare deployment files while building
[ -f "package-appengine.json" ] && cp package-appengine.json package.json

# Wait for build to complete
wait $BUILD_PID
[ ! -d "src/client/dist" ] && { echo -e "${RED}❌ Build failed${NC}"; exit 1; }

# Fast deployment
if [ "$ENV" = "dev" ]; then
  # Deploy to dev service (no promotion = faster)
  gcloud app deploy app-dev.yaml --no-promote --quiet --no-cache &
else
  # Deploy to production
  gcloud app deploy app.yaml --quiet --no-cache &
fi
DEPLOY_PID=$!

# Upload static files in parallel
BUCKET_NAME="${PROJECT_ID}-static"
gsutil -m -q cp -r src/client/dist/* gs://${BUCKET_NAME}/ &
STATIC_PID=$!

# Wait for both operations
wait $DEPLOY_PID
wait $STATIC_PID

# Cleanup
[ -f "package-appengine.json" ] && git checkout package.json 2>/dev/null || true

echo -e "${GREEN}✅ Deployment complete!${NC}"

if [ "$ENV" = "dev" ]; then
  DEV_URL="https://dev-dot-${PROJECT_ID}.appspot.com"
  echo -e "${YELLOW}🌐 Dev App: ${DEV_URL}${NC}"
  echo -e "${YELLOW}📦 Static: https://storage.googleapis.com/${BUCKET_NAME}/index.html${NC}"
  echo -e "${BLUE}💡 Custom domain: dev.queb.space${NC}"
  echo -e "${BLUE}   To set up: gcloud app domain-mappings create dev.queb.space --service=dev${NC}"
  
  # Auto-open dev environment (prefer custom domain)
  echo -e "${BLUE}🌐 Opening dev environment...${NC}"
  
  # Check if custom domain is available, otherwise use App Engine URL
  OPEN_URL="${DEV_URL}"
  if curl -s --connect-timeout 2 https://dev.queb.space > /dev/null 2>&1; then
    OPEN_URL="https://dev.queb.space"
    echo -e "${GREEN}✅ Using custom domain: dev.queb.space${NC}"
  fi
  
  if command -v open &> /dev/null; then
    open "${OPEN_URL}"
  elif command -v xdg-open &> /dev/null; then
    xdg-open "${OPEN_URL}"
  elif command -v start &> /dev/null; then
    start "${OPEN_URL}"
  else
    echo -e "${YELLOW}💡 Manually open: ${OPEN_URL}${NC}"
  fi
else
  PROD_URL="https://${PROJECT_ID}.appspot.com"
  echo -e "${YELLOW}🌐 Production App: ${PROD_URL}${NC}"
  echo -e "${YELLOW}📦 Static: https://storage.googleapis.com/${BUCKET_NAME}/index.html${NC}"
  echo -e "${BLUE}💡 Custom domain: queb.space${NC}"
  echo -e "${BLUE}   To set up: gcloud app domain-mappings create queb.space${NC}"
fi