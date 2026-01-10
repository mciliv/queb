# 🚀 Queb AI Setup - Environment-Agnostic Configuration

This guide ensures Queb's AI configuration works consistently across **any editor** (Cursor, Neovim, VSCode) and **any AI provider** (OpenAI, xAI).

## 📋 Quick Setup Checklist

### 1. Environment Variables
```bash
# Copy template
cp .env.example .env

# Edit .env with your API keys
# AI_PROVIDER=openai  # or 'xai'
# OPENAI_API_KEY=sk-...
# XAI_API_KEY=xai-...
```

### 2. Validate Configuration
```bash
npm run validate:ai
```

### 3. Start Development
```bash
npm start  # http://localhost:8080
```

## 🔒 Locked-In Configuration

### Environment Variables (`.env`)
```bash
# Required: Choose your AI provider
AI_PROVIDER=openai    # or 'xai'

# OpenAI Configuration (if using OpenAI)
OPENAI_API_KEY=sk-proj-...your-key-here...
OPENAI_MODEL=latest            # Auto-updates to newest GPT

# xAI Configuration (if using xAI)
XAI_API_KEY=xai-...your-key-here...
XAI_MODEL=latest               # Auto-updates to newest Grok
```

### Dependencies (Already Installed)
- `@ai-sdk/openai@^3.0.7` - Official OpenAI SDK
- `@ai-sdk/xai@^3.0.10` - Official xAI SDK
- `ai@^6.0.21` - Vercel AI SDK core
- `openai@^6.15.0` - OpenAI SDK (legacy support)

## 🧪 Validation

Run this command **anytime** to verify everything works:

```bash
npm run validate:ai
```

**Success Output:**
```
🎉 SUCCESS: AI Setup is LOCKED IN!
✅ Environment variables configured
✅ Dependencies installed
✅ AI Service working
✅ Provider switching functional
🚀 Ready to use with Cursor, Neovim, or any editor!
```

## 🔄 Provider Switching

### Via .env File (Primary Method)
```bash
# Edit .env file
AI_PROVIDER=xai    # Switch to xAI/Grok
AI_PROVIDER=openai # Switch to OpenAI/GPT

# Restart application
npm start
```

### Runtime Switching
```javascript
// Switch providers programmatically
const aiService = container.get('aiService');

aiService.switchProvider('xai');   // Switch to xAI (uses .env settings)
aiService.switchProvider('openai'); // Switch to OpenAI (uses .env settings)

// Same API for both providers
const result = await aiService.callAPI({
  messages: [{ role: 'user', content: 'Analyze this molecule...' }]
});
```

## 🛠️ Troubleshooting

### Validation Fails
```bash
npm run validate:ai  # Check specific errors
```

### Common Issues
- **Missing API keys**: Add them to `.env`
- **Wrong provider**: Set `AI_PROVIDER=openai` or `AI_PROVIDER=xai`
- **Network issues**: Check internet connection
- **Rate limits**: Wait and retry

### Manual Testing
```bash
# Test OpenAI
curl -X POST http://localhost:8080/api/structuralize \
  -H "Content-Type: application/json" \
  -d '{"text": "water molecule"}'

# Test xAI (switch provider first)
# Set AI_PROVIDER=xai in .env and restart
```

## 📁 File Structure (Locked-In)

```
/project/
├── .env                    # Your API keys (ignored by git)
├── .env.example           # Template (committed to git)
├── package.json           # Dependencies locked-in
├── scripts/
│   └── validate-ai-setup.js # Validation script
└── src/server/services/
    └── AIService.js       # Unified AI service
```

## 🎯 Editor Agnostic

This setup works identically in:
- **Cursor** - AI-powered IDE
- **Neovim** - Terminal-based editor
- **VSCode** - Microsoft's IDE
- **Any editor** - Just run `npm start`

## 🔑 API Keys

### OpenAI
1. Go to: https://platform.openai.com/api-keys
2. Create new secret key
3. Add to `.env`: `OPENAI_API_KEY=sk-proj-...`

### xAI
1. Go to: https://console.x.ai/
2. Create API key
3. Add to `.env`: `XAI_API_KEY=xai-...`

## 🚀 Production Deployment

The same `.env` structure works in production:

```bash
# Production environment variables
AI_PROVIDER=xai
OPENAI_API_KEY=...
XAI_API_KEY=...
NODE_ENV=production
```

## 💡 Pro Tips

- **Always run `npm run validate:ai`** before committing
- **Keep `.env` private** - Never commit API keys
- **Use "latest" models** for automatic updates
- **Test both providers** regularly
- **Document custom model changes** in commits

---

**Status**: ✅ **LOCKED IN** - Works across all editors and AI providers!