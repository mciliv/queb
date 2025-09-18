#!/bin/bash

# Check current environment status

echo "🔍 Environment Status:"
echo "===================="

# Check Node.js environment
if [ "$NODE_ENV" = "production" ]; then
    echo "🌐 Environment: PRODUCTION"
    echo "   🚨 Be careful with changes!"
else
    echo "💻 Environment: DEVELOPMENT (Local)"
    echo "   🛠️  Safe to experiment"
fi

# Check if server is running
if lsof -i :8080 > /dev/null 2>&1; then
    echo "🟢 Local Server: RUNNING (Port 8080)"
    echo "   📡 Access at: http://localhost:8080"
else
    echo "🔴 Local Server: NOT RUNNING"
    echo "   💡 Run 'd' to start local development"
fi

# Check if production domain is accessible
if curl -s --max-time 5 https://queb.space > /dev/null; then
    echo "🟢 Production: ACCESSIBLE"
    echo "   🌐 Live at: https://queb.space"
else
    echo "🟡 Production: NOT ACCESSIBLE"
    echo "   📝 May need deployment or DNS update"
fi

# Check git status
if [ -d ".git" ]; then
    echo ""
    echo "📋 Git Status:"
    git status --porcelain | head -5
    if [ $? -eq 0 ]; then
        echo "   ✅ Repository is clean"
    else
        echo "   📝 Uncommitted changes exist"
    fi
fi

echo ""
echo "🎯 Quick Actions:"
echo "   Local dev:  d"
echo "   Deploy:     ./scripts/deploy/deploy.sh"
echo "   Check logs: tail -f logs/server-*.log"


