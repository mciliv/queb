#!/bin/bash

# trust-cert.sh - Add development HTTPS certificate to system trust store

set -e

echo "🔒 Molecular Analysis App - Certificate Trust Setup"
echo "=================================================="

CERT_PATH="$(dirname "$0")/../../backend/api/certs/cert.pem"

if [ ! -f "$CERT_PATH" ]; then
    echo "❌ Certificate not found at: $CERT_PATH"
    echo "💡 Please start the development server first to generate certificates"
    echo "   Run: ./dev"
    exit 1
fi

echo "📄 Certificate found: $CERT_PATH"

# Detect OS and add certificate to trust store
case "$(uname -s)" in
    Darwin)
        echo "🍎 macOS detected - adding certificate to system keychain"
        sudo security add-trusted-cert -d -r trustRoot -k /Library/Keychains/System.keychain "$CERT_PATH"
        echo "✅ Certificate added to macOS system keychain"
        echo "🎉 You should no longer see 'trust website' prompts"
        ;;
    Linux)
        echo "🐧 Linux detected"
        if command -v update-ca-certificates >/dev/null; then
            # Ubuntu/Debian
            sudo cp "$CERT_PATH" /usr/local/share/ca-certificates/mol-dev.crt
            sudo update-ca-certificates
            echo "✅ Certificate added to Linux CA store"
        elif command -v update-ca-trust >/dev/null; then
            # CentOS/RHEL/Fedora
            sudo cp "$CERT_PATH" /etc/pki/ca-trust/source/anchors/mol-dev.crt
            sudo update-ca-trust
            echo "✅ Certificate added to Linux CA store"
        else
            echo "⚠️ Automatic certificate installation not supported"
            echo "💡 Please manually add this certificate to your browser:"
            echo "   Certificate: $CERT_PATH"
        fi
        ;;
    *)
        echo "⚠️ Unsupported OS: $(uname -s)"
        echo "💡 Please manually add this certificate to your browser:"
        echo "   Certificate: $CERT_PATH"
        ;;
esac

echo ""
echo "🌐 Alternative: Browser-specific setup"
echo "   Chrome: Visit https://localhost:3001 → Advanced → Proceed to localhost (unsafe)"
echo "   Firefox: Visit https://localhost:3001 → Advanced → Accept the Risk and Continue"
echo "   Safari: Visit https://localhost:3001 → Show Details → Visit this website"
echo ""
echo "🔄 Restart your browser after certificate installation" 