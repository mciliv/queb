#!/usr/bin/env node

const express = require('express');
const cors = require('cors');
const AIService = require('./src/server/services/agents/AIService');
const path = require('path');
const fs = require('fs');

// CLI argument parsing for configuration
const args = process.argv.slice(2);
const parsedArgs = {};

for (let i = 0; i < args.length; i++) {
  const arg = args[i];
  if (arg.startsWith('--')) {
    const [key, value] = arg.slice(2).split('=');
    parsedArgs[key] = value || true;
  }
}

// Extract parameters
const port = parsedArgs.port ? parseInt(parsedArgs.port) : (process.env.PORT || 8080);
const testMode = parsedArgs.test;
const devMode = parsedArgs.dev;

// Store configuration globally for endpoints to use
global.serverConfig = { testMode, devMode, port };

const app = express();

// Lazy initialization of AI Service (unified SDK abstraction)
// Trade-off evaluated: Keeping lazy init for dev/demo flexibility.
// Alternative (import at top) would validate config at startup but break server startup without API keys.
// Current approach allows server to run without AI features for demo purposes.
let aiService = null;
function getAIService() {
  if (!aiService) {
    aiService = new AIService();
  }
  return aiService;
}

// Middleware
app.use(cors());
app.use(express.json());
app.use(express.static('src/client/dist'));
app.use('/assets', express.static('src/client/dist/assets'));
app.use('/dist', express.static('src/client/dist'));

// Serve CSS file with correct name
app.get('/assets/style.css', (req, res) => {
  res.sendFile(path.join(__dirname, 'src', 'client', 'dist', 'bundle.css'));
});

// Simple chemical analysis endpoint
app.post('/api/analyze', async (req, res) => {
  try {
    const { text } = req.body;

    if (!text || typeof text !== 'string') {
      return res.status(400).json({ error: 'Text input required' });
    }

    // API key validation is now handled by AIService internally

    console.log(`Analyzing: "${text}"`);

    // #region agent log
    fetch('http://127.0.0.1:7243/ingest/f1225f0b-6c5b-477f-bc5d-1e74641debf9',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'server.js:36',message:'Analysis request received',data:{text,hasApiKey:!!process.env.OPENAI_API_KEY},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'D'})}).catch(()=>{});
    // #endregion

    const prompt = `Analyze what chemical compounds might be found in: ${text}

Return a JSON object with:
- object: the input text
- chemicals: array of objects with {name, smiles} for each chemical
- reason: brief explanation

Example for "coffee":
{
  "object": "coffee",
  "chemicals": [
    {"name": "Caffeine", "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"},
    {"name": "Chlorogenic acid", "smiles": "C1C(C(C(CC1(C(=O)O)O)OC(=O)C=CC2=CC(=C(C=C2)O)O)O)O"}
  ],
  "reason": "Coffee contains caffeine and various acids"
}`;

    const aiResponse = await getAIService().callAPI({
      messages: [{ role: 'user', content: prompt }],
      temperature: 0.1, // Low temperature for consistent results
      max_tokens: 1000
    });

    const response = aiResponse.content;

    // #region agent log
    fetch('http://127.0.0.1:7243/ingest/f1225f0b-6c5b-477f-bc5d-1e74641debf9',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'server.js:62',message:'AI response received',data:{responseLength:response?.length||0,startsWithBrace:response?.trim().startsWith('{')||false},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'C'})}).catch(()=>{});
    // #endregion

    // AIService automatically parses JSON responses - check if we got parsed object or raw content
    if (typeof aiResponse === 'object' && aiResponse !== null && !Array.isArray(aiResponse)) {
      // Already parsed JSON object from AIService
      const result = aiResponse;

      // #region agent log
      fetch('http://127.0.0.1:7243/ingest/f1225f0b-6c5b-477f-bc5d-1e74641debf9',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'server.js:65',message:'API response parsed',data:{hasChemicals:!!result.chemicals,chemicalCount:result.chemicals?.length||0,hasSmiles:result.chemicals?.some(c=>c.smiles)||false},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'C'})}).catch(()=>{});
      // #endregion

      res.json(result);
    } else {
      // Raw text response (AIService couldn't parse as JSON)
      // #region agent log
      fetch('http://127.0.0.1:7243/ingest/f1225f0b-6c5b-477f-bc5d-1e74641debf9',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'server.js:69',message:'Non-JSON response',data:{responseLength:response?.length||0},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'C'})}).catch(()=>{});
      // #endregion

      // If AI returns non-JSON, wrap it
      res.json({
        object: text,
        chemicals: [],
        reason: response || 'Analysis completed',
        raw_response: response
      });
    }

  } catch (error) {
    console.error('Analysis error:', error);
    res.status(500).json({
      error: 'Analysis failed',
      details: error.message
    });
  }
});

// Health check
app.get('/api/health', (req, res) => {
  res.json({
    status: 'ok',
    timestamp: new Date().toISOString(),
    config: global.serverConfig
  });
});

// Test endpoints
app.get('/api/test', (req, res) => {
  console.log('🧪 Running tests via endpoint...');
  const { spawn } = require('child_process');

  const testProcess = spawn('node', ['test-runner.js'], {
    cwd: process.cwd()
  });

  let output = '';
  let errorOutput = '';

  testProcess.stdout.on('data', (data) => {
    output += data.toString();
  });

  testProcess.stderr.on('data', (data) => {
    errorOutput += data.toString();
  });

  testProcess.on('close', (code) => {
    res.json({
      success: code === 0,
      code,
      output,
      error: errorOutput
    });
  });
});

app.get('/api/test/:file', (req, res) => {
  const testFile = req.params.file;
  console.log(`🧪 Running specific test: ${testFile}`);

  try {
    // Try to run the test file directly
    require(`./${testFile}`);
    res.json({ success: true, message: `Test ${testFile} completed` });
  } catch (error) {
    res.status(500).json({
      success: false,
      error: error.message,
      stack: error.stack
    });
  }
});

// Serve the main HTML page (React client) with basic template replacement
app.get('/', (req, res) => {
  try {
    const indexPath = path.join(__dirname, 'src', 'client', 'dist', 'index.html');
    let htmlContent = fs.readFileSync(indexPath, 'utf8');

    // Basic template replacement
    const cacheBust = Date.now();
    htmlContent = htmlContent
      .split('{{CACHE_BUST}}')
      .join(cacheBust);

    res.type('html').send(htmlContent);
  } catch (error) {
    console.error('Failed to serve index.html:', error);
    res.status(500).send('Server error');
  }
});

// Start server
app.listen(PORT, () => {
  console.log(`🚀 Chemical Analyzer running on http://localhost:${PORT}`);
  console.log(`📝 Enter text like "coffee", "water", "aspirin" to analyze chemical content`);
});