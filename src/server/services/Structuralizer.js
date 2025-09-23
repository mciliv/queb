const { createClient } = require("../ai/openai/client");
const aiConfig = require("../ai/config");
const logger = require("./file-logger");
const promptEngine = require('../../core/PromptEngine');
const errorHandler = require('../../core/ErrorHandler');
const MolecularProcessor = require("./molecular-processor");
const { resolveName, getPropertiesByCID } = require("./name-resolver");

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

  }

  // Pattern according to https://platform.openai.com/docs/overview
  // as of 2025-09-21
  async callOpenAI(requestParams) {
    if (!this.isOpenAIAvailable) throw new Error("AI service unavailable");
    
    const activeModel = this.resolvedModelName || this.model;
    
    // Modern OpenAI API call
    const response = await this.client.chat.completions.create({
      model: activeModel,
      ...requestParams
    });
    
    this.resolvedModelName = activeModel;
    return response;
  }

  // Structuralize multimodal: accepts { object?, imageBase64?, x?, y? } and returns { object, chemicals, recommendedBox?, reason? }
  async structuralize(payload) {
    const inputObject = typeof payload?.object === 'string' ? payload.object.trim() : '';
    const imageBase64 = typeof payload?.imageBase64 === 'string' ? payload.imageBase64 : null;
    const x = typeof payload?.x === 'number' ? payload.x : null;
    const y = typeof payload?.y === 'number' ? payload.y : null;

    let objectText = inputObject;
    let recommendedBox = null;
    let reason = null;

    if (!objectText && imageBase64) {
      const byImage = await this.structuralizeImage(imageBase64, null, x, y, null, null, null).catch(() => null);
      if (byImage) {
        objectText = byImage.object || '';
        recommendedBox = byImage.recommendedBox || null;
        reason = byImage.reason || null;
      }
    }

    const byText = await this.structuralizeText(objectText || '');
    return {
      object: byText.object || objectText || (imageBase64 ? 'image' : ''),
      chemicals: Array.isArray(byText.chemicals) ? byText.chemicals : [],
      recommendedBox: recommendedBox || null,
      reason: reason || byText.reason || null
    };
  }

  // Primary: text → molecules (names→structures) using a single structuralization prompt
  async structuralizeText(object) {
    logger.info(`🔍 Starting text analysis for: "${object}"`);
    
    if (!this.isOpenAIAvailable) {
      logger.error("❌ AI service unavailable for structuralization");
      throw new Error("AI service unavailable for structuralization");
    }
    
    const prompt = this.testConfig.prompt && this.testConfig.prompt !== 'custom'
      ? this.testConfig.prompt
      : promptEngine.generateChemicalPrompt(object || '', { includeReason: true });
    
    logger.info(`📝 Generated prompt for analysis`);
    
    let parsed = { object, chemicals: [] }; // Initialize with fallback
    try {
      // Modern OpenAI API call for structured analysis
      const response = await this.callOpenAI({
        messages: [{ 
          role: 'user', 
          content: prompt 
        }],
        response_format: { type: 'json_object' }
      });
      
      logger.info(`✅ Received AI response`);
      const content = response.choices[0].message.content;
      logger.info(`📄 Raw AI response: ${content?.substring(0, 200)}...`);
      
      // Use PromptEngine for JSON parsing and repair
      parsed = promptEngine.repairJSON(content);
      if (parsed) {
        logger.info(`✅ Successfully parsed JSON response`);
      } else {
        logger.warn(`⚠️ Failed to parse JSON response, using fallback`);
        parsed = { object, chemicals: [] };
      }
    } catch (aiError) {
      const handled = errorHandler.handleAIError(aiError, { method: 'structuralizeText', object });
      logger.error(`❌ AI call failed:`, handled);
      throw new Error(handled.message);
    }
    const list = Array.isArray(parsed?.chemicals) ? parsed.chemicals : [];
    // Programmatic fallback: fill missing SMILES (and CID) using resolvers
    const enriched = [];
    for (const item of list) {
      const name = item?.name || '';
      let smiles = item?.smiles || null;
      let cid = item?.cid ?? null;
      if (!smiles && name) {
        try {
          if (cid) {
            const props = await getPropertiesByCID(cid).catch(() => null);
            smiles = (props && props.smiles) ? props.smiles : null;
          }
          if (!smiles) {
            const res = await resolveName(name).catch(() => null);
            cid = (res && res.cid) ? res.cid : cid;
            smiles = (res && res.smiles) ? res.smiles : smiles;
          }
        } catch (resolveError) {
          logger.warn(`⚠️ Failed to resolve name "${name}":`, resolveError.message);
        }
      }
      enriched.push({ name, smiles: smiles || null, cid: cid ?? null });
    }

    // Attempt SDF generation for identified molecules by name or smiles (after enrichment)
    const resolved = [];
    for (const item of enriched) {
      const name = item?.name || '';
      const smiles = item?.smiles || null;
      let sdfPath = null;
      try {
        if (smiles) {
          sdfPath = await this.molecularProcessor.generateSDF(smiles, false).catch(() => null);
        }
        if (!sdfPath && name) {
          const byName = await this.molecularProcessor.generateSDFByName(name, false).catch(() => null);
          if (byName && byName.sdfPath) sdfPath = byName.sdfPath;
        }
      } catch (sdfError) {
        logger.warn(`⚠️ Failed to generate SDF for "${name}" (${smiles}):`, sdfError.message);
      }
      resolved.push({ name, smiles: smiles || null, sdfPath: sdfPath || null, status: sdfPath ? 'ok' : (smiles ? 'smiles_only' : 'lookup_required') });
    }
    
    logger.info(`🎯 Text analysis complete: Found ${resolved.length} chemicals`);
    resolved.forEach((chem, i) => {
      logger.info(`   ${i+1}. ${chem.name} - Status: ${chem.status}`);
    });
    
    return { object: parsed.object || object, chemicals: resolved };
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
    try {
      if (!content || typeof content !== "string") return { object: "Unknown object", chemicals: [] };
      try { return JSON.parse(content); } catch (_) {}
      const match = content.match(/\{[\s\S]*\}/); if (match) { try { return JSON.parse(match[0]); } catch (_) {} }
      return { object: "Unknown object", chemicals: [] };
    } catch (_) { return { object: "Unknown object", chemicals: [] }; }
  }
}

module.exports = Structuralizer;


