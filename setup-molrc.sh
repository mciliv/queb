#!/bin/bash

# Setup script for .molrc system
# Run once to configure the molecular project

echo "🧬 Setting up Molecular Project Configuration System..."

# Check Python environment
echo "🐍 Checking Python environment..."
if command -v python3 >/dev/null 2>&1; then
    PYTHON_VERSION=$(python3 --version 2>&1)
    echo "✅ Found: $PYTHON_VERSION"
    
    # Check for required Python packages
    echo "📦 Checking Python packages..."
    python3 -c "import rdkit; print('✅ RDKit available')" 2>/dev/null || echo "⚠️  RDKit not found - chemistry tests may fail"
    python3 -c "import pytest; print('✅ pytest available')" 2>/dev/null || echo "⚠️  pytest not found - Python tests may fail"
else
    echo "⚠️  Python3 not found - Python functionality will not work"
fi

# Check Node environment
echo "📦 Checking Node.js environment..."
if command -v node >/dev/null 2>&1; then
    NODE_VERSION=$(node --version)
    echo "✅ Found Node.js: $NODE_VERSION"
    
    if command -v npm >/dev/null 2>&1; then
        NPM_VERSION=$(npm --version)
        echo "✅ Found npm: $NPM_VERSION"
    else
        echo "⚠️  npm not found"
    fi
else
    echo "❌ Node.js not found - project will not work"
fi

# Create alias in shell profile
SHELL_PROFILE=""
if [ -f ~/.zshrc ]; then
    SHELL_PROFILE="$HOME/.zshrc"
elif [ -f ~/.bashrc ]; then
    SHELL_PROFILE="$HOME/.bashrc"
elif [ -f ~/.bash_profile ]; then
    SHELL_PROFILE="$HOME/.bash_profile"
fi

if [ -n "$SHELL_PROFILE" ]; then
    # Check if alias already exists
    if ! grep -q "alias mol=" "$SHELL_PROFILE"; then
        echo "" >> "$SHELL_PROFILE"
        echo "# Molecular Project Commands" >> "$SHELL_PROFILE"
        echo "alias mol='./molrc'" >> "$SHELL_PROFILE"
        echo "✅ Added 'mol' alias to $SHELL_PROFILE"
        echo "   Run 'source $SHELL_PROFILE' or restart terminal to use"
    else
        echo "✅ 'mol' alias already exists"
    fi
else
    echo "⚠️  Could not find shell profile to add alias"
fi

echo ""
echo "🎯 Usage:"
echo "  mol              # Show project info and available scripts"
echo "  mol list         # List all scripts"
echo "  mol dev          # 🚀 DO EVERYTHING: clean, build, test, browser, server"
echo "  mol build        # Build frontend only"
echo "  mol test         # Run tests only"
echo "  mol edit         # Edit .molrc configuration"
echo ""
echo "📝 To add new scripts:"
echo "  1. Edit .molrc file"
echo "  2. Add to [scripts] section"
echo "  3. Run 'mol list' to verify"
echo ""
echo "✅ Setup complete! This is now the definitive script registry."
