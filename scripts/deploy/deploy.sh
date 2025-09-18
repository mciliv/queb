#!/bin/bash

# Deploy to queb.space
# This script handles production deployment

set -e

echo "🚀 Deploying to queb.space..."
echo "📝 Committing local changes..."

# Add all changes
git add .

# Commit with timestamp
git commit -m "Deploy $(date +'%Y-%m-%d %H:%M:%S')" || echo "No changes to commit"

# Push to remote
echo "📤 Pushing to remote repository..."
git push origin main

# Deploy to production (you'll need to configure this based on your hosting)
echo "🌐 Deploying to queb.space..."
echo "   Configure your deployment method:"
echo "   - Vercel: vercel --prod"
echo "   - Netlify: netlify deploy --prod"
echo "   - Railway: railway up"
echo "   - GCP: gcloud app deploy"
echo ""
echo "✅ Code pushed! Configure your hosting provider to deploy from the repository."
echo "🌐 Your app will be available at: https://queb.space"


