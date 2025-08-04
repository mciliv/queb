#!/bin/bash

# Simple Project Initialization
# ============================

# Load toolkit
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../toolkit/utils.sh"

# Setup logging
setup_logging "$(basename "$0")"

PROJECT_ROOT="$(dirname "$(dirname "$0")")"
PROJECT_NAME="$(basename "$PROJECT_ROOT")"

log_info "Initializing project: $PROJECT_NAME"

# Check if we're in a git repository
if ! git rev-parse --git-dir > /dev/null 2>&1; then
    log_error "Not in a git repository. Please run this from your project root."
    exit 1
fi

# Check if project is already initialized
if [ -f "scripts/config.sh" ]; then
    log_warn "Project appears to be already initialized"
    echo "Options:"
    echo "1. Reinitialize (will overwrite existing config)"
    echo "2. Skip initialization"
    read -p "Choose option (1-2): " choice
    
    case $choice in
        1)
            log_info "Reinitializing project..."
            rm -f scripts/config.sh
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
1. Clone the repository
2. Initialize the project:
   \`\`\`bash
   ./scripts/init-simple.sh
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

### Deployment
\`\`\`bash
./toolkit/deploy.sh
\`\`\`

### Domain Check
\`\`\`bash
source toolkit/utils.sh
source scripts/config.sh
check_domain "\$DOMAIN_NAME"
\`\`\`

## Scripts

This project uses a simple toolkit for deployment operations.
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
    "test": "jest",
    "deploy": "./toolkit/deploy.sh"
  },
  "keywords": [],
  "author": "",
  "license": "ISC",
  "dependencies": {
    "express": "^4.18.2",
    "cors": "^2.8.5"
  },
  "devDependencies": {
    "jest": "^29.0.0",
    "nodemon": "^3.0.0"
  }
}
EOF
    log_info "Created package.json"
fi

# Update the config with project name
if [ -f "scripts/config.sh" ]; then
    log_info "Updating configuration with project name..."
    sed -i.bak "s/PROJECT_NAME=\"your-project\"/PROJECT_NAME=\"$PROJECT_NAME\"/" scripts/config.sh
    rm scripts/config.sh.bak 2>/dev/null || true
fi

# Create example usage script
EXAMPLE_SCRIPT="scripts/example-usage.sh"
if [ ! -f "$EXAMPLE_SCRIPT" ]; then
    cat > "$EXAMPLE_SCRIPT" << 'EOF'
#!/bin/bash

# Example Usage Script
# ===================
# Shows how to use the simplified toolkit

# Load toolkit
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../toolkit/utils.sh"
source "$SCRIPT_DIR/config.sh"

# Setup logging
setup_logging "$(basename "$0")"

log_info "Running example script for $PROJECT_NAME"

# Use toolkit functions
check_auth
check_domain "$DOMAIN_NAME"
check_python "$PYTHON_VERSION"

log_info "Example script completed successfully"
EOF
    chmod +x "$EXAMPLE_SCRIPT"
    log_info "Created example script: $EXAMPLE_SCRIPT"
fi

# Commit all changes
log_info "Committing project initialization..."
git add .
git commit -m "Initialize project with simplified toolkit"

log_info "✅ Project initialization completed!"
echo ""
echo "Next steps:"
echo "1. Update scripts/config.sh with your project-specific values"
echo "2. Set up your environment variables"
echo "3. Install dependencies: npm install"
echo "4. Test the setup: ./scripts/example-usage.sh"
echo "5. Deploy: ./toolkit/deploy.sh" 