const { createClient } = require("../ai/openai/client");
const fs = require('fs');
const path = require('path');
const aiConfig = require("../ai/config");
const logger = require("./file-logger");
const promptEngine = require('../../core/PromptEngine');
const errorHandler = require('../../core/ErrorHandler');
const MolecularProcessor = require("./molecular-processor");
<<<<<<< Updated upstream
const { resolveName, getPropertiesByCID } = require("./name-resolver");
const { convertNamesToSmiles } = require("../prompts/name-to-smiles");
const wikidata = require('./wikidata-resolver');

class Structuralizer {
  constructor(apiKey = null, testConfig = null) {
    // Single provider/model
    this.client = createClient();
    const effectiveApiKey = apiKey || aiConfig.apiKey;
    this.isOpenAIAvailable = !!(this.client && effectiveApiKey);
    this.model = aiConfig.model;
    this.resolvedModelName = null;
    this.chemicalInstructions = null;
    this.molecularProcessor = new MolecularProcessor();

    // Test configuration support
    this.testConfig = testConfig || {};
    if (this.testConfig.model) {
      this.model = this.testConfig.model;
    }
=======
const { resolveName } = require("./name-resolver");
>>>>>>> Stashed changes

/**
 * Main function: Analyzes any input modality (text, image, or both) and returns structured chemical data
 * @param {Object} payload - Input with { object?, imageBase64?, x?, y?, lookupMode?, testConfig? }
 * @returns {Promise<Object>} - { object, chemicals, recommendedBox?, reason? }
 */
async function chemicals(payload) {
  // Memory check at start
  const memStart = process.memoryUsage();
  if (memStart.heapUsed > 400 * 1024 * 1024) { // 400MB threshold
    throw new Error("Server memory limit exceeded - please try again later");
  }

  const inputObject = typeof payload?.object === 'string' ? payload.object.trim() : '';
  const imageBase64 = typeof payload?.imageBase64 === 'string' ? payload.imageBase64 : null;
  const x = typeof payload?.x === 'number' ? payload.x : null;
  const y = typeof payload?.y === 'number' ? payload.y : null;
  const lookupMode = payload?.lookupMode || 'database';
  const testConfig = payload?.testConfig || {};

  // Initialize AI client
  const client = createClient();
  const apiKey = aiConfig.apiKey;
  const isAIAvailable = !!(client && apiKey);
  const model = testConfig.model || aiConfig.model;
  const molecularProcessor = new MolecularProcessor();

  if (!isAIAvailable) {
    logger.error("❌ AI service unavailable - client or API key missing");
    throw new Error("AI service unavailable - check OpenAI API key configuration");
  }

  // Helper: Call OpenAI API
  const callAI = async (requestParams) => {
    const timeoutMs = Math.min(10000, Number(process.env.AI_TIMEOUT_MS || 10000));
    const aiPromise = client.chat.completions.create({ model, ...requestParams });
    const timeoutPromise = new Promise((_, reject) => {
      setTimeout(() => reject(new Error('timeout')), timeoutMs);
    });
    return await Promise.race([aiPromise, timeoutPromise]);
  };

  let objectText = inputObject;
  let recommendedBox = null;
  let reason = null;

  // Step 1: If image provided and no text, extract object name from image
  if (!objectText && imageBase64) {
    logger.info(`🖼️ Analyzing image to detect object`);
    try {
      let imageData = imageBase64;
      if (/^https?:\/\//i.test(imageBase64)) {
        const fetch = require('node-fetch');
        const resp = await fetch(imageBase64);
        if (!resp.ok) throw new Error(`Failed to fetch image URL: ${resp.status}`);
        const buf = await resp.buffer();
        imageData = buf.toString('base64');
      }

      // Load detection prompt base text and append click coordinates if provided
      const detectionBasePath = path.join(__dirname, '../../core/prompts/object_detection.txt');
      let detectionPrompt = '';
      try { detectionPrompt = fs.readFileSync(detectionBasePath, 'utf8').trim(); } catch (_) { detectionPrompt = ''; }
      if (Number.isFinite(x) && Number.isFinite(y)) {
        detectionPrompt += `\n\nUser clicked at coordinates (${x}, ${y}). Focus analysis around this area.`;
      }

      const messages = [{ 
        role: 'user', 
        content: [
          { type: 'text', text: detectionPrompt },
          { type: 'image_url', image_url: { url: `data:image/jpeg;base64,${imageData}`, detail: 'high' } }
        ] 
      }];

      const response = await callAI({ messages, max_tokens: 500, response_format: { type: 'json_object' } });
      const content = response.choices[0].message.content;
      let parsed;
      try { parsed = JSON.parse(content); } catch (_) { parsed = { object: 'Unknown object' }; }
      
      objectText = parsed.object || 'Unknown object';
      reason = typeof parsed?.reason === 'string' ? parsed.reason : null;

      // Extract bounding box if provided
      const box = parsed?.recommendedBox || parsed?.box || parsed?.recommendedCrop || {};
      const bx = Number(box.x), by = Number(box.y), bw = Number(box.width), bh = Number(box.height);
      if ([bx, by, bw, bh].every(n => Number.isFinite(n)) && bw > 0 && bh > 0) {
        recommendedBox = { x: Math.round(bx), y: Math.round(by), width: Math.round(bw), height: Math.round(bh) };
      }

      logger.info(`✅ Image analysis complete: "${objectText}"`);
    } catch (imageError) {
      logger.warn(`⚠️ Image analysis failed:`, imageError.message);
      objectText = 'Material in image';
    }
  }

  // Step 2: Analyze text (from input or extracted from image) to identify chemicals
  logger.info(`🔍 Starting chemical analysis for: "${objectText}"`);
  
  const prompt = (testConfig.prompt && testConfig.prompt !== 'custom')
    ? testConfig.prompt
    : (() => {
        // Load minimal chemical analysis prompt and append object + optional reason suffix
        const chemBasePath = path.join(__dirname, '../../core/prompts/chemical_analysis.txt');
        const reasonPath = path.join(__dirname, '../../core/prompts/reason_suffix.txt');
        let base = '';
        let reasonSuffix = '';
        try { base = fs.readFileSync(chemBasePath, 'utf8').trim(); } catch (_) { base = ''; }
        try { reasonSuffix = fs.readFileSync(reasonPath, 'utf8').trim(); } catch (_) { reasonSuffix = ''; }
        return [
          base,
          `Object: ${objectText || ''}`,
          reasonSuffix
        ].filter(Boolean).join('\n\n');
      })();
  
  let parsed = { object: objectText, chemicals: [] };
  try {
    const response = await callAI({
      messages: [{ role: 'user', content: prompt }],
      response_format: { type: 'json_object' }
    });
    
    const content = response.choices[0].message.content;
    logger.info(`📄 Raw AI response: ${content?.substring(0, 200)}...`);
    
    parsed = promptEngine.repairJSON(content);
    if (!parsed) {
      logger.warn(`⚠️ Failed to parse JSON response, using fallback`);
      parsed = { object: objectText, chemicals: [] };
    } else {
      logger.info(`✅ Successfully parsed JSON response`);
    }
  } catch (aiError) {
    logger.warn(`⚠️ AI call failed, using fallback:`, aiError?.message || String(aiError));
    const handled = errorHandler.handleAIError(aiError, { method: 'chemicals', object: objectText });
    reason = handled.message || 'AI unavailable';
    parsed = { object: objectText, chemicals: [] };
  }

  const list = Array.isArray(parsed?.chemicals) ? parsed.chemicals : [];
  logger.info(`🔧 Using lookup mode: ${lookupMode}`);
  
  // Step 3: Enrich SMILES based on lookup mode
  const enriched = [];
  for (const item of list) {
    const name = item?.name || '';
    let smiles = item?.smiles || null;
    
    if (lookupMode === 'database') {
      // Database-first: Always try database lookup
      if (name) {
        try {
          const res = await resolveName(name).catch(() => null);
          smiles = (res && res.smiles) ? res.smiles : (item?.smiles || null);
          logger.debug(`Database lookup for "${name}": ${smiles ? 'SUCCESS' : 'FAILED'}`);
        } catch (resolveError) {
          logger.warn(`⚠️ Database lookup failed for "${name}":`, resolveError.message);
          smiles = item?.smiles || null;
        }
      }
    } else {
      // AI-first: Use AI SMILES, database fills gaps
      if (!smiles && name) {
        try {
          const res = await resolveName(name).catch(() => null);
          smiles = (res && res.smiles) ? res.smiles : smiles;
          logger.debug(`Database fallback for "${name}": ${smiles ? 'SUCCESS' : 'FAILED'}`);
        } catch (resolveError) {
          logger.warn(`⚠️ Database fallback failed for "${name}":`, resolveError.message);
        }
      }
    }
    
    enriched.push({ name, smiles: smiles || null });
  }

  // Step 4: Generate SDF files for molecules
  const resolved = [];
  for (const item of enriched) {
    const name = item?.name || '';
    const smiles = item?.smiles || null;
    let sdfPath = null;
    
<<<<<<< Updated upstream
    return { object: parsed.object || object, chemicals: resolved };
  }

  // Stage 1: return only object specification (no molecules)
  async specifyObject(payload) {
    const inputObject = typeof payload?.object === 'string' ? payload.object.trim() : '';
    const imageBase64 = typeof payload?.imageBase64 === 'string' ? payload.imageBase64 : null;
    const x = typeof payload?.x === 'number' ? payload.x : null;
    const y = typeof payload?.y === 'number' ? payload.y : null;

    let objectText = inputObject;
    let recommendedBox = null;
    let reason = null;

    // If object not provided but image present, detect object via vision
    if (!objectText && imageBase64) {
      const byImage = await this.structuralizeImage(imageBase64, null, x, y, null, null, null).catch(() => null);
      if (byImage) {
        objectText = byImage.object || '';
        recommendedBox = byImage.recommendedBox || null;
        reason = byImage.reason || null;
      }
    }

    const objectType = this.detectObjectType(objectText || '');
    return { object: objectText || (imageBase64 ? 'image' : ''), objectType, recommendedBox, reason };
  }

  // Stage 2: resolve components using DB-first, then LLM-assisted names listing
  async resolveComponents(object) {
    const objectText = typeof object === 'string' ? object.trim() : '';
    if (!objectText) return { object: '', molecules: [] };

    // DB-first (Object-level): Try Wikidata for real-world object constituents
    try {
      const parts = await wikidata.findObjectConstituents(objectText);
      if (Array.isArray(parts) && parts.length > 0) {
        // Map to names then resolve via DB for SMILES/CID
        const names = parts.map(p => ({ name: p.label })).filter(p => p.name);
        if (names.length > 0) {
          const resolved = await convertNamesToSmiles({ object: objectText, molecules: names }, this.client).catch(() => ({ molecules: [] }));
          if (resolved && Array.isArray(resolved.molecules) && resolved.molecules.length > 0) {
            return { object: objectText, molecules: resolved.molecules };
          }
        }
      }
    } catch (_) {}

    // DB-first (Direct compound): if the object itself is a specific compound, resolve directly
    try {
      const direct = await resolveName(objectText);
      if (direct && (direct.smiles || direct.cid)) {
        // Heuristic: if PubChem returns a direct match with SMILES or CID and the title closely matches the object
        const title = (direct.title || '').toLowerCase();
        const normalizedObj = objectText.toLowerCase();
        const isDirectMatch = title === normalizedObj || title.includes(normalizedObj) || normalizedObj.includes(title);
        if (isDirectMatch && (direct.smiles || direct.cid)) {
          const molecules = [{ name: direct.title || objectText, cid: direct.cid || null, smiles: direct.smiles || null, status: (direct.smiles ? 'ok' : 'lookup_required') }];
          return { object: objectText, molecules };
        }
      }
    } catch (_) {}

    // LLM names listing (no SMILES), then resolve via PubChem/DB, then optional LLM fallback for SMILES
    if (!this.isOpenAIAvailable) {
      // If AI unavailable, return empty and let caller handle
      return { object: objectText, molecules: [] };
    }

    // Ask LLM to propose specific molecule names only
    const namePrompt = require('../../core/PromptEngine').generateNamePrompt(objectText);
    let namesOnly = [];
    try {
      const response = await this.callOpenAI({
        messages: [{ role: 'user', content: namePrompt }],
        response_format: { type: 'json_object' }
      });
      const content = response.choices[0].message.content;
      const parsed = require('../../core/PromptEngine').repairJSON(content) || {};
      const mols = Array.isArray(parsed.molecules) ? parsed.molecules : [];
      namesOnly = mols.map(m => ({ name: typeof m?.name === 'string' ? m.name.trim() : '' })).filter(m => m.name);
    } catch (_) {
      namesOnly = [];
    }

    if (namesOnly.length === 0) {
      // As a last resort, return empty set and let client decide next steps
      return { object: objectText, molecules: [] };
    }

    // Resolve names to CIDs/SMILES using DB-first then LLM fallback (convertNamesToSmiles already does this)
    const resolved = await convertNamesToSmiles({ object: objectText, molecules: namesOnly }, this.client).catch(() => ({ object: objectText, molecules: [] }));
    return { object: resolved.object || objectText, molecules: Array.isArray(resolved.molecules) ? resolved.molecules : [] };
  }

  // Secondary: object detection from image → object text (+ optional bounding box)
  async structuralizeImage(imageBase64, croppedImageBase64 = null, x = null, y = null, cropMiddleX = null, cropMiddleY = null, cropSize = null) {
    if (!imageBase64 || imageBase64.trim().length === 0) throw new Error("Empty image data");
    if (!this.isOpenAIAvailable) throw new Error("AI service unavailable");
    if (/^https?:\/\//i.test(imageBase64)) {
      const fetch = require('node-fetch');
      const resp = await fetch(imageBase64);
      if (!resp.ok) throw new Error(`Failed to fetch image URL: ${resp.status}`);
      const buf = await resp.buffer();
      imageBase64 = buf.toString('base64');
    }
    const detectionText = promptEngine.generateDetectionPrompt({ x, y });
    const messages = [{ role: 'user', content: [
      { type: 'text', text: detectionText },
      { type: 'image_url', image_url: { url: `data:image/jpeg;base64,${imageBase64}`, detail: 'high' } }
    ] }];
    if (croppedImageBase64) {
      messages[0].content.push({ type: 'text', text: `Here is a cropped view near the focus area.` });
      messages[0].content.push({ type: 'image_url', image_url: { url: `data:image/jpeg;base64,${croppedImageBase64}`, detail: 'high' } });
    }
    const response = await this.callOpenAI({ messages, max_tokens: 500, response_format: { type: 'json_object' } });
    const content = response.choices[0].message.content;
    let parsed; try { parsed = JSON.parse(content); } catch (_) { parsed = { object: 'Unknown object' }; }
    const box = parsed?.recommendedBox || parsed?.box || parsed?.recommendedCrop || {};
    const bx = Number(box.x); const by = Number(box.y); const bw = Number(box.width); const bh = Number(box.height);
    let recommendedBox = null;
    if ([bx, by, bw, bh].every(n => Number.isFinite(n)) && bw > 0 && bh > 0) {
      recommendedBox = { x: Math.round(bx), y: Math.round(by), width: Math.round(bw), height: Math.round(bh) };
    }
    const reason = typeof parsed?.reason === 'string' ? parsed.reason : null;
    return { object: parsed.object || 'Unknown object', recommendedBox, reason };
  }

  // Removed legacy extract/convert helpers (simplified to two-step flow)

  detectObjectType(object) {
    const objectLower = object.toLowerCase();
    if (objectLower.includes('wine') || objectLower.includes('beer') || objectLower.includes('coffee') || objectLower.includes('drink') || objectLower.includes('beverage')) return 'beverage';
    if (objectLower.includes('fruit') || objectLower.includes('food') || objectLower.includes('apple') || objectLower.includes('orange') || objectLower.includes('vegetable')) return 'food';
    if (objectLower.includes('plastic') || objectLower.includes('metal') || objectLower.includes('wood') || objectLower.includes('stone')) return 'material';
    return 'general';
  }

  parseAIResponse(content) {
=======
>>>>>>> Stashed changes
    try {
      const withTimeout = (p, ms) => new Promise((resolve) => {
        let done = false;
        const t = setTimeout(() => { if (!done) resolve(null); }, ms);
        p.then((v) => { if (!done) { done = true; clearTimeout(t); resolve(v); } })
         .catch(() => { if (!done) { done = true; clearTimeout(t); resolve(null); } });
      });

      if (smiles) {
        const res = await withTimeout(molecularProcessor.generateSDF(smiles, false), 3000);
        if (res && res.sdfPath) sdfPath = res.sdfPath;
      }
      if (!sdfPath && name) {
        const byName = await withTimeout(molecularProcessor.generateSDFByName(name, false), 3000);
        if (byName && byName.sdfPath) sdfPath = byName.sdfPath;
      }
    } catch (sdfError) {
      logger.warn(`⚠️ Failed to generate SDF for "${name}" (${smiles}):`, sdfError.message);
    }
    
    resolved.push({ 
      name, 
      smiles: smiles || null, 
      sdfPath: sdfPath || null, 
      status: sdfPath ? 'ok' : (smiles ? 'smiles_only' : 'lookup_required') 
    });
  }
  
  logger.info(`🎯 Chemical analysis complete: Found ${resolved.length} chemicals`);
  resolved.forEach((chem, i) => {
    logger.info(`   ${i+1}. ${chem.name} - Status: ${chem.status}`);
  });
  
  // Memory cleanup and final check
  const memEnd = process.memoryUsage();
  const memUsed = memEnd.heapUsed - memStart.heapUsed;
  logger.info(`Memory usage: ${Math.round(memUsed/1024/1024)}MB increase`);
  
  // Force garbage collection if available
  if (global.gc) {
    global.gc();
  }
  
  return {
    object: parsed.object || objectText || (imageBase64 ? 'Material in image' : ''),
    chemicals: resolved,
    recommendedBox: recommendedBox || null,
    reason: reason || parsed.reason || null
  };
}

// Backward compatibility: Export as class for existing code
// DEPRECATED: This class will be removed in v2.0.0
class Structuralizer {
  constructor(apiKey = null, testConfig = null) {
    console.warn(
      '[DEPRECATED] Structuralizer class is deprecated and will be removed in v2.0.0.\n' +
      'Migration: Use the chemicals() function instead:\n' +
      '  Before: const s = new Structuralizer(apiKey); await s.structuralizeText("aspirin");\n' +
      '  After:  const { chemicals } = require("./Structuralizer"); await chemicals({ object: "aspirin" });'
    );
    this.testConfig = testConfig || {};
  }
  
  async structuralize(payload) {
    logger.warn('[DEPRECATED] structuralize() → Use chemicals() directly');
    return chemicals({ ...payload, testConfig: this.testConfig });
  }
  
  async structuralizeText(object, lookupMode = 'database') {
    logger.warn('[DEPRECATED] structuralizeText() → Use chemicals({ object, lookupMode })');
    return chemicals({ object, lookupMode, testConfig: this.testConfig });
  }
  
  async structuralizeImage(imageBase64, croppedImageBase64 = null, x = null, y = null) {
    logger.warn('[DEPRECATED] structuralizeImage() → Use chemicals({ imageBase64, x, y })');
    return chemicals({ imageBase64, x, y, testConfig: this.testConfig });
  }
}

// Export the new API as primary
module.exports = { 
  chemicals,
  Structuralizer  // Deprecated - for backward compatibility only
};

// For CommonJS default import compatibility
module.exports.default = chemicals;


