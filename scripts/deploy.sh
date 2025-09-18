#!/bin/bash

# Molecular Analysis App Deployment Script
# Deploys to multiple targets based on environment

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}🚀 Starting deployment...${NC}"

# Check if we're in the right directory
if [ ! -f "package.json" ]; then
    echo -e "${RED}❌ Error: package.json not found. Are you in the project root?${NC}"
    exit 1
fi

# Parse deployment target
DEPLOY_TARGET=${1:-"local"}
echo -e "${BLUE}📍 Deployment target: ${DEPLOY_TARGET}${NC}"

# Build the project
echo -e "${BLUE}🔨 Building project...${NC}"
npm run build

# Check if build was successful
if [ ! -d "src/client/dist" ]; then
    echo -e "${RED}❌ Error: Build failed - dist directory not found${NC}"
    exit 1
fi

echo -e "${GREEN}✅ Build completed successfully${NC}"

# Deployment logic based on target
case $DEPLOY_TARGET in
    "local")
        echo -e "${BLUE}📁 Local deployment - copying files to local server directory${NC}"
        LOCAL_DEPLOY_DIR="${HOME}/deployed-apps/molecular-analysis"
        mkdir -p "$LOCAL_DEPLOY_DIR"
        
        # Copy built files
        cp -r src/client/dist/* "$LOCAL_DEPLOY_DIR/"
        cp -r public/* "$LOCAL_DEPLOY_DIR/" 2>/dev/null || true
        
        echo -e "${GREEN}✅ Local deployment completed${NC}"
        echo -e "${YELLOW}📄 Files deployed to: $LOCAL_DEPLOY_DIR${NC}"
        ;;
        
    "staging"|"stage")
        echo -e "${BLUE}🎭 Staging deployment${NC}"
        echo -e "${YELLOW}💡 Configure your staging server details in this script${NC}"
        echo -e "${YELLOW}   Example: rsync -av src/client/dist/ user@staging-server:/var/www/molecular-analysis/${NC}"
        ;;
        
    "production"|"prod")
        echo -e "${BLUE}🌟 Production deployment${NC}"
        
        # Safety check for production
        echo -e "${YELLOW}⚠️  This will deploy to PRODUCTION. Continue? (y/N)${NC}"
        read -r confirmation
        if [[ $confirmation != [yY] ]]; then
            echo -e "${YELLOW}❌ Production deployment cancelled${NC}"
            exit 0
        fi
        
        echo -e "${YELLOW}💡 Configure your production deployment in this script${NC}"
        echo -e "${YELLOW}   Examples:${NC}"
        echo -e "${YELLOW}   - Cloud Functions: gcloud functions deploy molecular-analysis --source .${NC}"
        echo -e "${YELLOW}   - Vercel: vercel --prod${NC}"
        echo -e "${YELLOW}   - Netlify: netlify deploy --prod --dir=src/client/dist${NC}"
        ;;
        
    "docker")
        echo -e "${BLUE}🐳 Docker deployment${NC}"
        
        # Build Docker image
        if [ -f "Dockerfile" ]; then
            docker build -t molecular-analysis:latest .
            echo -e "${GREEN}✅ Docker image built${NC}"
            echo -e "${YELLOW}💡 To run: docker run -p 8080:8080 molecular-analysis:latest${NC}"
        else
            echo -e "${YELLOW}💡 Create a Dockerfile to enable Docker deployment${NC}"
        fi
        ;;
        
    *)
        echo -e "${RED}❌ Unknown deployment target: $DEPLOY_TARGET${NC}"
        echo -e "${YELLOW}💡 Available targets: local, staging, production, docker${NC}"
        echo -e "${YELLOW}   Usage: ship [target]${NC}"
        echo -e "${YELLOW}   Example: ship production${NC}"
        exit 1
        ;;
esac

# Post-deployment checks
echo -e "${BLUE}🔍 Running post-deployment checks...${NC}"

# Check if files exist
if [ "$DEPLOY_TARGET" = "local" ] && [ -f "$LOCAL_DEPLOY_DIR/bundle.js" ]; then
    echo -e "${GREEN}✅ Bundle files deployed successfully${NC}"
fi

# Display next steps
echo -e "${GREEN}🎉 Deployment completed successfully!${NC}"
echo ""
echo -e "${BLUE}📋 Next steps:${NC}"
case $DEPLOY_TARGET in
    "local")
        echo -e "${YELLOW}   • Files are in: $LOCAL_DEPLOY_DIR${NC}"
        echo -e "${YELLOW}   • Set up a web server to serve these files${NC}"
        ;;
    "staging"|"production")
        echo -e "${YELLOW}   • Configure your server deployment commands in this script${NC}"
        echo -e "${YELLOW}   • Set up environment variables on your server${NC}"
        ;;
    "docker")
        echo -e "${YELLOW}   • Push to registry: docker push your-registry/molecular-analysis:latest${NC}"
        echo -e "${YELLOW}   • Deploy to your container platform${NC}"
        ;;
esac

echo -e "${BLUE}📄 To customize deployment: edit $0${NC}"