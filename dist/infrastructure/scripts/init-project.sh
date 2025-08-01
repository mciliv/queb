#!/bin/bash

# Initialize Project Script
# ========================
# Sets up a new project with mol-toolkit and basic structure

# Fallback logging (toolkit doesn't exist yet)
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}
log_info() { log "INFO: $1"; }
log_error() { log "ERROR: $1"; }
log_warn() { log "WARN: $1"; }

# Configuration
TOOLKIT_REPO="https://github.com/mciliv/mol-toolkit.git"
PROJECT_ROOT="$(dirname "$(dirname "$0")")"
PROJECT_NAME="$(basename "$PROJECT_ROOT")"

log_info "Initializing project: $PROJECT_NAME"

# Check if we're in a git repository
if ! git rev-parse --git-dir > /dev/null 2>&1; then
    log_error "Not in a git repository. Please run this from your project root."
    exit 1
fi

# Check if project is already initialized
if [ -d ".mol-toolkit" ] || [ -f "scripts/config.sh" ]; then
    log_warn "Project appears to be already initialized"
    echo "Options:"
    echo "1. Reinitialize (will overwrite existing config)"
    echo "2. Skip initialization"
    read -p "Choose option (1-2): " choice
    
    case $choice in
        1)
            log_info "Reinitializing project..."
            rm -rf .mol-toolkit scripts/config.sh scripts/project-specific
            ;;
        2)
            log_info "Skipping initialization"
            exit 0
            ;;
        *)
            log_error "Invalid choice"
            exit 1
            ;;
    esac
fi

# Run the setup script
log_info "Setting up mol-toolkit..."
if [ -f "scripts/setup-toolkit.sh" ]; then
    ./scripts/setup-toolkit.sh
else
    log_error "setup-toolkit.sh not found. Please run this from a project with the setup script."
    exit 1
fi

# Create basic project structure
log_info "Creating project structure..."

# Create logs directory
mkdir -p logs

# Create .gitignore entries for logs
if [ ! -f ".gitignore" ]; then
    cat > .gitignore << 'EOF'
# Logs
logs/
*.log

# Environment variables
.env
.env.local

# Dependencies
node_modules/
__pycache__/
*.pyc

# Build outputs
dist/
build/

# IDE
.vscode/
.idea/

# OS
.DS_Store
Thumbs.db
EOF
    log_info "Created .gitignore"
fi

# Create basic README if it doesn't exist
if [ ! -f "README.md" ]; then
    cat > README.md << EOF
# $PROJECT_NAME

## Overview
[Add your project description here]

## Setup

### Prerequisites
- Node.js
- Python 3.12+
- Google Cloud CLI
- Git

### Installation
1. Clone the repository:
   \`\`\`bash
   git clone --recurse-submodules <your-repo-url>
   \`\`\`

2. Initialize the project:
   \`\`\`bash
   ./scripts/init-project.sh
   \`\`\`

3. Update configuration:
   \`\`\`bash
   # Edit scripts/config.sh with your project-specific values
   \`\`\`

4. Set up environment variables:
   \`\`\`bash
   export OPENAI_API_KEY="your-api-key"
   export GOOGLE_CLOUD_PROJECT="your-project-id"
   \`\`\`

## Usage

### Domain Management
\`\`\`bash
./scripts/project-specific/check-domain-status.sh
\`\`\`

### Python Environment
\`\`\`bash
source scripts/project-specific/helper.sh
\`\`\`

### Deployment
\`\`\`bash
npm run deploy-to-production
\`\`\`

## Scripts

See \`scripts/README.md\` for detailed documentation of available scripts.

## Toolkit

This project uses the mol-toolkit for common operations. See \`.mol-toolkit/README.md\` for toolkit documentation.
EOF
    log_info "Created README.md"
fi

# Create package.json if it doesn't exist
if [ ! -f "package.json" ]; then
    cat > package.json << 'EOF'
{
  "name": "your-project-name",
  "version": "1.0.0",
  "description": "Your project description",
  "main": "server.js",
  "scripts": {
    "start": "node server.js",
    "dev": "nodemon server.js",
    "test": "npm run test:unit && npm run test:integration",
    "test:unit": "jest --testPathPattern=unit.test.js --verbose --silent",
    "test:integration": "jest --testPathPattern=integration.test.js --verbose",
    "deploy-to-production": "node deploy-to-production",
    "ship-to-production-with-git": "node ship-to-production-with-git"
  },
  "keywords": [],
  "author": "",
  "license": "ISC",
  "dependencies": {
    "express": "^4.18.2",
    "cors": "^2.8.5",
    "openai": "^4.0.0"
  },
  "devDependencies": {
    "jest": "^29.0.0",
    "nodemon": "^3.0.0",
    "supertest": "^6.3.0"
  }
}
EOF
    log_info "Created package.json"
fi

# Update the config with project name
if [ -f "scripts/config.sh" ]; then
    log_info "Updating configuration with project name..."
    sed -i.bak "s/PROJECT_NAME=\"your-project-name\"/PROJECT_NAME=\"$PROJECT_NAME\"/" scripts/config.sh
    rm scripts/config.sh.bak 2>/dev/null || true
fi

# Create example project-specific script
EXAMPLE_SCRIPT="scripts/project-specific/example.sh"
if [ ! -f "$EXAMPLE_SCRIPT" ]; then
    cat > "$EXAMPLE_SCRIPT" << 'EOF'
#!/bin/bash

# Example Project-Specific Script
# ===============================
# This shows how to use the toolkit in your project

# Load toolkit utilities
TOOLKIT_DIR="$(dirname "$(dirname "$(dirname "${BASH_SOURCE[0]}")")")/.mol-toolkit"
source "${TOOLKIT_DIR}/core/logging-utils.sh"
source "${TOOLKIT_DIR}/core/gcloud-utils.sh"

# Load project configuration
source "$(dirname "$0")/../config.sh"

# Setup error handling and logging
set_error_handling "$(basename "$0")"

# Main script logic
log_info "Running example script for $PROJECT_NAME"
log_info "Project: $PROJECT_NAME"
log_info "Region: $REGION"

# Use toolkit functions
if check_gcloud_auth; then
    log_info "Google Cloud is authenticated"
else
    log_error "Google Cloud not authenticated"
    exit 1
fi

log_info "Example script completed successfully"
EOF
    chmod +x "$EXAMPLE_SCRIPT"
    log_info "Created example script: $EXAMPLE_SCRIPT"
fi

# Commit all changes
log_info "Committing project initialization..."
git add .
git commit -m "Initialize project with mol-toolkit and basic structure"

log_info "✅ Project initialization completed!"
echo ""
echo "Next steps:"
echo "1. Update scripts/config.sh with your project-specific values"
echo "2. Set up your environment variables (OPENAI_API_KEY, etc.)"
echo "3. Install dependencies: npm install"
echo "4. Test the setup: ./scripts/project-specific/example.sh"
echo ""
echo "Documentation:"
echo "- Project scripts: scripts/README.md"
echo "- Toolkit: .mol-toolkit/README.md" 