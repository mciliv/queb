#!/bin/bash

# Quick Environment Check
source "$(dirname "$0")/utils.sh"

echo "🔍 Environment Status"
echo "===================="

# Environment
[ "$NODE_ENV" = "production" ] && echo "🚨 PRODUCTION" || echo "💻 DEVELOPMENT"

# Local server
if lsof -i :8080 > /dev/null 2>&1; then
    echo "🟢 Server: RUNNING (http://localhost:8080)"
else
    echo "🔴 Server: NOT RUNNING (run 'd' to start)"
fi

# Production check
if curl -s --max-time 5 https://queb.space > /dev/null; then
    echo "🟢 Production: ACCESSIBLE (https://queb.space)"
else
    echo "🟡 Production: NOT ACCESSIBLE"
fi

# Git status
if [ -d ".git" ]; then
    changes=$(git status --porcelain | wc -l)
    [ "$changes" -eq 0 ] && echo "✅ Git: Clean" || echo "📝 Git: $changes changes"
fi

echo ""
echo "Quick commands: d (dev) | d deploy (deploy) | d clean (cleanup)"


