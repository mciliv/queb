#!/bin/bash

# Setup script for .molrc system
# Run once to configure the molecular project

echo "🧬 Setting up Molecular Project Configuration System..."

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
