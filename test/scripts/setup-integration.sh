#!/bin/bash

# Integration Test Setup Script
# Ensures all dependencies are ready for integration testing

set -e

echo "🔧 Setting up integration test environment..."

# Check if required tools are installed
command -v node >/dev/null 2>&1 || { echo "❌ Node.js is required but not installed"; exit 1; }
command -v npm >/dev/null 2>&1 || { echo "❌ npm is required but not installed"; exit 1; }

# Install dependencies if needed
if [ ! -d "node_modules" ]; then
    echo "📦 Installing dependencies..."
    npm install
fi

# Check for Python dependencies if chemistry tests are needed
if [ -f "chemistry/processors/sdf.py" ]; then
    echo "🐍 Checking Python dependencies..."
    python3 -c "import rdkit" 2>/dev/null || echo "⚠️ RDKit not found - chemistry tests may fail"
fi

# Ensure test directories exist
mkdir -p test/results
mkdir -p data/sdf_files
mkdir -p logs

# Kill any existing processes on test ports
echo "🧹 Cleaning up existing processes..."
lsof -ti:8080 | xargs kill -9 2>/dev/null || true
lsof -ti:3001 | xargs kill -9 2>/dev/null || true

# Setup test database (if needed)
if [ -f "infrastructure/scripts/setup-database.sh" ]; then
    echo "🗄️ Setting up test database..."
    NODE_ENV=test ./infrastructure/scripts/setup-database.sh
fi

echo "✅ Integration test environment ready!"
echo ""
echo "Available commands:"
echo "  npm run test:unit                    # Run unit tests only"
echo "  npm run test:integration            # Run integration tests"
echo "  npm run test:watch                  # Watch files and run tests"
echo "  node test/config/test-runner.js --all # Run all tests in order"
echo ""