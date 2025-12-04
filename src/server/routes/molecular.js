// Molecular analysis and PubChem API endpoints
const express = require('express');
const router = express.Router();
const path = require('path');
const fs = require('fs');
const crypto = require('crypto');

const setupMolecularRoutes = (chemicals, molecularProcessor, resolveName) => {
  /*
   * Function Diagram:
   * 
   * setupMolecularRoutes(chemicals, molecularProcessor, resolveName)
   * │
   * ├─── Utility Functions
   * │    ├─ ensureSdfDir()           → Creates/returns SDF directory path
   * │    ├─ sanitizeName(raw)        → Cleans/hashes names for filenames
   * │    ├─ saveSdf(filename, text)  → Writes SDF file, returns URL path
   * │    ├─ fetchText(url)           → Generic HTTP text fetcher
   * │    └─ fetchPubchemSdf({...})   → Fetches SDF from PubChem API
   * │
   * ├─── Route Handlers
   * │    │
   * │    ├─ POST /structuralize
   * │    │   ├─ Validates request body
   * │    │   ├─ Checks memory limits (500MB threshold)
   * │    │   ├─ Calls → chemicals(req.body)
   * │    │   ├─ Tracks duration & memory usage
   * │    │   └─ Returns → processed chemical data
   * │    │
   * │    ├─ POST /name-to-smiles
   * │    │   ├─ Takes array of chemical names
   * │    │   ├─ Calls → resolveName(name) for each
   * │    │   └─ Returns → [{name, smiles, status}]
   * │    │
   * │    ├─ POST /name-to-sdf
   * │    │   ├─ Takes single chemical name
   * │    │   ├─ Calls → molecularProcessor.generateSDFByName(name, overwrite)
   * │    │   └─ Returns → {name, sdfPath, status}
   * │    │
   * │    ├─ POST /pubchem/sdf
   * │    │   ├─ Takes name or smiles + record_type
   * │    │   ├─ Calls → fetchPubchemSdf({...})
   * │    │   ├─ Calls → saveSdf() to persist
   * │    │   └─ Returns → {sdfPath, status, source}
   * │    │
   * │    └─ POST /generate-sdfs
   * │        ├─ Validates with SdfGenerationSchema
   * │        ├─ Takes array of SMILES strings
   * │        ├─ Calls → molecularProcessor.processSmiles(smiles, overwrite)
   * │        └─ Returns → {sdfPaths, errors, skipped, message}
   * │
   * └─── Returns → router (Express Router instance)
   * 
   * Data Flow:
   *   External Request → Route Handler (validates input)
   *                   → Injected Dependency (chemicals/molecularProcessor/resolveName)
   *                   → Utility Functions (if file operations needed)
   *                   → Response to Client
   */
  
  // Utility functions for SDF handling
  const ensureSdfDir = () => {
    const dir = path.join(
      __dirname,
      "..",
      "..",
      process.env.NODE_ENV === "test" ? "tests/sdf_files" : "public/sdf_files"
    );
    if (!fs.existsSync(dir)) fs.mkdirSync(dir, { recursive: true });
    return dir;
  };

  const sanitizeName = (raw) => {
    try {
      const base = String(raw || '').trim();
      if (!base) return 'unknown';
      const cleaned = base.replace(/[^a-zA-Z0-9]+/g, '_').replace(/^_+|_+$/g, '');
      if (cleaned.length <= 64) return cleaned;
      const hash = crypto.createHash('md5').update(base).digest('hex').slice(0, 12);
      return `${cleaned.slice(0, 32)}_${hash}`;
    } catch (_) {
      return 'unknown';
    }
  };

  const saveSdf = (filename, text) => {
    const dir = ensureSdfDir();
    const file = path.join(dir, filename);
    fs.writeFileSync(file, text, 'utf8');
    return `/sdf_files/${filename}`;
  };

  const fetchText = async (url) => {
    try {
      if (typeof fetch !== 'undefined') {
        const resp = await fetch(url);
        if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
        return await resp.text();
      }
    } catch (_) {}
    try {
      const fetchLib = require('node-fetch');
      const resp = await fetchLib(url);
      if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
      return await resp.text();
    } catch (_) {
      throw new Error('No fetch available');
    }
  };

  const fetchPubchemSdf = async ({ name, smiles, record_type = '3d' }) => {
    let selected = null;
    if (name) {
      selected = { type: 'name', value: name };
    } else if (smiles) {
      selected = { type: 'smiles', value: smiles };
    }

    if (!selected) {
      throw new Error('Provide name or smiles');
    }

    const enc = encodeURIComponent(selected.value);
    const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/${selected.type}/${enc}/SDF?record_type=${record_type}`;

    const resp = await fetch(url);
    if (!resp.ok) {
      throw new Error(`PubChem fetch failed: HTTP ${resp.status}`);
    }

    return resp.text();
  };

  /**
   * 
   * This is the primary API endpoint that powers Queb
   * It accepts various input types and returns identified chemical compounds.
   * 
   * Request flow:
   * 1. Validate request payload
   * 2. Check server memory (prevent OOM)
   * 3. Call chemicals() service for analysis
   * 4. Track performance metrics
   * 5. Return results with chemical data and SDF paths
   * 
   * Request body:
   * - object: Text description (e.g., "coffee")
   * - imageBase64: Base64 image data
   * - x, y: Click coordinates for image analysis
   * - lookupMode: "database" or "ai"
   * 
   * Response:
   * - object: Identified object name
   * - chemicals: Array of {name, sdfPath, status}
   * - recommendedBox: Suggested crop region for images
   */
  router.post("/structuralize", async (req, res) => {
    // Set timeout to prevent hanging requests (30 seconds)
    const timeout = setTimeout(() => {
      if (!res.headersSent) {
        res.status(408).json({ error: "Request timeout - processing took too long" });
      }
    }, 30000);

    try {
      // Step 1: Validate request payload
      if (!req.body || (typeof req.body !== 'object')) {
        clearTimeout(timeout);
        return res.status(400).json({ error: "Invalid payload" });
      }

      // Step 2: Check memory usage to prevent server crashes
      const memBefore = process.memoryUsage();
      if (memBefore.heapUsed > 500 * 1024 * 1024) {  // 500MB limit
        clearTimeout(timeout);
        return res.status(507).json({ error: "Server memory limit exceeded" });
      }

      // Step 3: Perform molecular analysis
      const startedAt = Date.now();
      const out = await chemicals(req.body);  // Main analysis function
      const durationMs = Date.now() - startedAt;

      // Step 4: Track performance metrics for monitoring
      const memAfter = process.memoryUsage();
      const memIncrease = memAfter.heapUsed - memBefore.heapUsed;

      // Log analysis results for debugging and monitoring
      try {
        const obj = out?.object || req.body?.object || '-';
        const list = Array.isArray(out?.chemicals) ? out.chemicals : (Array.isArray(out?.molecules) ? out.molecules : []);
        const count = list.length;
        const sample = list.slice(0, 3);  // Log first 3 chemicals as sample
        console.log(`[structuralize] object="${obj}" count=${count} duration=${durationMs}ms mem=${Math.round(memIncrease/1024/1024)}MB sample=`, sample);
      } catch (_) {}  // Logging should never break the response

      // Step 5: Send successful response
      clearTimeout(timeout);
      res.json(out);
    } catch (error) {
      clearTimeout(timeout);
      console.error('[structuralize] failed:', error?.message || error);
      res.status(500).json({ error: `Structuralization failed: ${error.message}` });
    }
  });

  /**
   * POST /name-to-smiles - Convert chemical names to SMILES notation
   * 
   * Batch endpoint for resolving multiple chemical names to SMILES.
   * SMILES (Simplified Molecular Input Line Entry System) is needed
   * for 3D structure generation.
   * 
   * Request: { items: [{name: "caffeine"}, {name: "aspirin"}] }
   * Response: { molecules: [{name, smiles, status}] }
   */
  router.post("/name-to-smiles", async (req, res) => {
    try {
      const { items } = req.body || {};
      if (!Array.isArray(items)) {
        return res.status(400).json({ error: "items array is required" });
      }
      
      // Process each chemical name
      const results = [];
      for (const it of items) {
        let name = it?.name || '';
        let smiles = null;
        
        // Attempt to resolve name to SMILES via PubChem
        try {
          if (!smiles) {
            const res = await resolveName(name).catch(() => null);
            if (res?.smiles) {
              smiles = res.smiles;
            }
          }
        } catch (_) {}  // Continue with null if resolution fails
        
        results.push({ 
          name, 
          smiles, 
          status: smiles ? 'ok' : 'lookup_required' 
        });
      }
      
      res.json({ molecules: results });
    } catch (error) {
      res.status(500).json({ error: error.message });
    }
  });

  router.post("/name-to-sdf", async (req, res) => {
    try {
      const { name, overwrite = false } = req.body || {};
      if (typeof name !== 'string' || name.trim().length === 0) {
        return res.status(400).json({ error: "name is required" });
      }

      const byName = await molecularProcessor.generateSDFByName(name, overwrite);
      if (!byName || !byName.sdfPath) {
        return res.json({ name, sdfPath: null, status: 'lookup_required' });
      }

      res.json({ name: byName.name || name, sdfPath: byName.sdfPath, status: 'ok' });
    } catch (error) {
      res.status(500).json({ error: error.message });
    }
  });

  router.post('/pubchem/sdf', async (req, res) => {
    try {
      const { name, smiles, record_type = '3d' } = req.body || {};
      if ((!name || name.trim().length === 0) && (!smiles || smiles.trim().length === 0)) {
        return res.status(400).json({ error: 'Provide name or smiles' });
      }
      const text = await fetchPubchemSdf({ name, smiles, record_type });
      const safe = sanitizeName(`${name || smiles}_${record_type}`);
      const sdfPath = saveSdf(`${safe}.sdf`, text);
      res.json({ sdfPath, status: 'ok', source: 'pubchem' });
    } catch (error) {
      if (error.message.includes('404') || error.message.includes('NOT_FOUND')) {
        return res.status(404).json({ error: 'not found from provided identifiers' });
      }
      res.status(500).json({ error: error.message || 'fetch failed' });
    }
  });

  router.post("/generate-sdfs", async (req, res) => {
    try {
      const { SdfGenerationSchema } = require("../schemas/schemas");
      const validation = SdfGenerationSchema.safeParse(req.body);
      if (!validation.success) {
        return res.status(400).json({
          error: "Invalid input data",
          details: validation.error.issues,
        });
      }

      const { smiles, overwrite = false } = req.body;

      if (!smiles || !Array.isArray(smiles)) {
        return res.status(400).json({ error: "smiles array is required" });
      }

      const result = await molecularProcessor.processSmiles(smiles, overwrite);

      res.json({
        sdfPaths: result.sdfPaths,
        errors: result.errors,
        skipped: result.skipped,
        message: "Files generated",
      });
    } catch (error) {
      res.status(500).json({ error: error.message });
    }
  });

  // Text-based molecular analysis endpoints
  router.post('/object-molecules', async (req, res) => {
    try {
      const { object } = req.body || {};
      if (!object || typeof object !== 'string' || object.trim().length === 0) {
        return res.status(400).json({ error: 'object string is required' });
      }

      // Use the chemicals function (Structuralizer) to analyze the object
      const result = await chemicals({ object: object.trim() });
      
      res.json({
        output: {
          object: object.trim(),
          chemicals: result.chemicals || result.molecules || []
        }
      });
    } catch (error) {
      console.error('[object-molecules] failed:', error?.message || error);
      res.status(500).json({ error: `Analysis failed: ${error.message}` });
    }
  });

  router.post('/structures-from-text', async (req, res) => {
    try {
      const { object } = req.body || {};
      if (!object || typeof object !== 'string' || object.trim().length === 0) {
        return res.status(400).json({ error: 'object string is required' });
      }

      // Use the chemicals function (Structuralizer) to analyze the object
      const result = await chemicals({ object: object.trim() });
      
      // Return direct format (no output wrapper)
      res.json({
        chemicals: result.chemicals || result.molecules || []
      });
    } catch (error) {
      console.error('[structures-from-text] failed:', error?.message || error);
      res.status(500).json({ error: `Analysis failed: ${error.message}` });
    }
  });

  return router;
};

module.exports = setupMolecularRoutes;