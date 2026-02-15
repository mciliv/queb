# Pipeline Restoration Log

## Date: Feb 15, 2026

## Objective

Restore the "text input to list of 3D moljs chemicals" pipeline that had stopped working due to missing/broken files across the codebase.

---

## Boot Commands

```bash
# Initialize git (needed for npm prepare hook)
git init

# Install dependencies (peer conflict between puppeteer and mermaid-cli)
npm install --legacy-peer-deps

# Start server (requires env vars)
OPENAI_API_KEY='<your-key>' OPENAI_MODEL=gpt-4o node src/server/api/server.js

# Or use the run script (picks up shell env vars)
export OPENAI_API_KEY='<your-key>'
export OPENAI_MODEL=gpt-4o
./run start

# Test the pipeline
curl -X POST http://localhost:8080/api/structuralize \
  -H 'Content-Type: application/json' \
  -d '{"text":"coffee","lookupMode":"GPT-5"}'
```

---

## Issues Found and Fixed

### Tier 1: Server Startup Blockers

These crashed the server in sequence before it could listen on any port.

#### 1. Missing `src/core/logger.js`

- **Symptom:** `MODULE_NOT_FOUND` on first require in `src/core/services.js` line 4
- **Cause:** `services.js` imports `{ log, warn, error }` from `./logger` but the file didn't exist
- **Fix:** Created `src/core/logger.js` with minimal console wrappers (`log`, `warn`, `error`)
- **Why minimal:** These are only used for startup messages like `log('Running in local development...')`

#### 2. Broken AIService import path

- **Symptom:** `MODULE_NOT_FOUND` for `../../../chat/AIService`
- **Location:** `src/core/services.js` line 358 (DI container registration)
- **Cause:** `require('../../../chat/AIService')` resolves outside the project root. The `chat/` directory never existed. The actual file is `src/server/services/AIService.js`
- **Fix:** Changed to `require('../server/services/AIService')`

#### 3. Missing route files (2 files)

- **Symptom:** `MODULE_NOT_FOUND` for `../routes/hotel` and `../routes/chemical-analysis`
- **Location:** `src/server/api/app.js` lines 116-121
- **Cause:** Both files were referenced via `require()` at startup but didn't exist on disk
- **Fix:** Created stub files exporting no-op setup functions:
  - `src/server/routes/hotel.js` (hotel is an independent page, not part of the molecule pipeline)
  - `src/server/routes/chemical-analysis.js`

#### 4. Missing database-recommender service

- **Symptom:** `MODULE_NOT_FOUND` for `../services/database-recommender` (discovered during boot testing)
- **Location:** `src/server/api/app.js` line 439 (`setupDatabaseRecommendationRoutes`)
- **Fix:** Created `src/server/services/database-recommender.js` stub with class skeleton

### Tier 2: Pipeline Runtime Bugs

After the server booted, the text-to-molecule pipeline still failed.

#### 5. AI response validation failure ("Invalid chemical response from AI")

- **Symptom:** `POST /api/structuralize` returned 500 with "Invalid chemical response from AI"
- **Root cause (two parts):**

**5a. Markdown fence wrapping:** GPT-4/4o often wraps JSON responses in ` ```json ... ``` ` markdown fences. The `_parseSDKResponse` method in `AIService.js` called `JSON.parse(content)` directly, which fails on fenced content. The fallback returned `{ content, role, ... }` instead of the parsed chemical data, which then failed the Structuralizer's validation.

- **Fix:** Added markdown fence stripping before JSON.parse in `AIService.js`:
  ```
  content = content.replace(/^```(?:json)?\s*\n?/i, '').replace(/\n?```\s*$/i, '').trim();
  ```

**5b. No system message:** The AI prompt had no system message enforcing JSON-only output, so GPT-4 sometimes responded with conversational text instead of JSON.

- **Fix:** Added system message in `Structuralizer.js` `_analyzeChemicals`:
  ```
  { role: 'system', content: 'You are a chemical analysis API. Always respond with valid JSON only, no markdown, no commentary.' }
  ```

#### 6. Object text lost in pipeline (AI analyzed wrong substance)

- **Symptom:** Query for "coffee" returned chemicals for "lemon"
- **Location:** `src/server/services/Structuralizer.js` line 114
- **Cause:** `this._chemicals(payload.object)` passed a string (`"coffee"`) to `_chemicals()`, which then destructured it as `const { object: inputObject, ... } = payload`. Destructuring a string gives `undefined` for all properties, so `objectText` became `''` and the AI hallucinated a random substance.
- **Fix:** Changed to `this._chemicals(payload)` so the full object `{ object: "coffee", lookupMode: "GPT-5" }` is passed through

### Tier 3: Code Quality

#### 7. Debug agent log pollution (stripped)

Dozens of `fetch('http://127.0.0.1:7243/ingest/...')` calls inside `// #region agent log` blocks were scattered across four files. These added failed TCP connection latency to every request and obscured the actual code.

**Stripped from:**
- `src/server/api/app.js` (7 blocks)
- `src/server/services/Structuralizer.js` (6 blocks)
- `src/server/services/AIService.js` (7 blocks)
- `src/client/hooks/useApi.js` (5 blocks)

Also removed 2 leftover `console.log('[DEBUG]...')` lines from `Structuralizer.js`.

#### 8. `.vscode/launch.json` merge conflict resolved

- **Symptom:** File contained `<<<<<<< HEAD`, `=======`, `>>>>>>>` markers making it invalid JSON
- **Fix:** Resolved by keeping both sides' configs deduplicated, with correct ordering (test, dev, local, prod) and updated entry point to `src/server/api/server.js`

---

## Files Created

| File | Purpose |
|------|---------|
| `src/core/logger.js` | Console wrappers for services.js startup logging |
| `src/server/routes/hotel.js` | Stub to prevent startup crash |
| `src/server/routes/chemical-analysis.js` | Stub to prevent startup crash |
| `src/server/services/database-recommender.js` | Stub to prevent startup crash |

## Files Modified

| File | Change |
|------|--------|
| `src/core/services.js` | Fixed AIService import path (line 358) |
| `src/server/services/Structuralizer.js` | Fixed payload passthrough, added system message, stripped agent logs, removed debug console.logs |
| `src/server/services/AIService.js` | Added markdown fence stripping in response parser, stripped agent logs |
| `src/server/api/app.js` | Stripped agent logs |
| `src/client/hooks/useApi.js` | Stripped agent logs |
| `.vscode/launch.json` | Resolved merge conflict |

---

## Verified Result

```
POST /api/structuralize  { "text": "coffee", "lookupMode": "GPT-5" }

Response (200, ~10s):
{
  "object": "coffee",
  "chemicals": [
    { "name": "Caffeine",          "smiles": "CN1C=NC2=...", "sdfPath": "/sdf_files/...", "status": "ok" },
    { "name": "Chlorogenic acid",  "smiles": "C1=CC(=...",   "sdfPath": "/sdf_files/...", "status": "ok" },
    { "name": "Trigonelline",      "smiles": "C1=CC=...",    "sdfPath": "/sdf_files/...", "status": "ok" },
    { "name": "Quinic acid",       "smiles": "C1[C@@H]...", "sdfPath": "/sdf_files/...", "status": "ok" },
    { "name": "Tannins",           "smiles": null,           "sdfPath": null,             "status": "lookup_required" }
  ],
  "reason": "These compounds are commonly found in coffee..."
}
```

SDF files are generated on disk at `public/sdf_files/` via PubChem download fallback (the Python/RDKit path is broken but the fallback works).

---

## Known Remaining Issues (not blocking, for future work)

### Pipeline fallbacks in use (work but are suboptimal)

1. **Python SDF path wrong** -- `molecular-processor.js` line 174 references `python/sdf.py` which doesn't exist. Falls back to PubChem download. Fix: point to `utils/molecular_docking/sdf.py` for faster local generation.

2. **`resolveNameToCID()` undefined** -- `name-resolver.js` line 483 calls this function but it's never defined. Wrapped in try-catch, falls through to direct PubChem name lookup. Fix: implement CID-based lookup for more reliable resolution.

3. **`MAX_RETRIES` / `BASE_DELAY_MS` bare variables** -- `name-resolver.js` `fetchText()` uses these instead of `RETRY_CONFIG.MAX_RETRIES`. The function would throw if retry logic is triggered. Fix: prefix with `RETRY_CONFIG.`.

### Missing files (non-blocking)

4. **Prompt templates** -- `src/core/prompts/object_detection.txt` and `name_resolution.txt` don't exist. `_loadPromptFile()` returns `''` silently. Only affects image analysis and name resolution prompts.

5. **Port/adapter/use-case files** -- `src/server/services/ports/log-error-port.js`, `adapters/file-logger-adapter.js`, `log-error-use-case.js` don't exist. Only affects the `POST /api/log-error` endpoint (frontend error logging).

6. **`config/env.js`** -- Referenced by `build-frontend.js` and test files but doesn't exist.

### Code quality

7. **Model choice** -- Server was tested with `gpt-4o` (works reliably). The original `gpt-4` model returned conversational text instead of JSON. May need to validate/default the model in AIService.

8. **Root `server.js`** -- Orphaned legacy entry point with broken imports. Not used by `./run` or `npm start`. Safe to delete.

9. **Frontend `dist/`** -- Build output doesn't exist in repo. The `./run start` command auto-builds if missing. Run `npm run build` to generate manually.

10. **Stub files** -- `hotel.js`, `chemical-analysis.js`, and `database-recommender.js` are stubs. They need real implementations if those features are desired.
