// Inline full implementation (renamed from AtomPredictor)
let OpenAIClient = null;
try { OpenAIClient = require("openai").OpenAI; } catch (_) { OpenAIClient = null; }
const { imageToSmiles_instructions, textToNames_prompt } = require("../prompts");
const MolecularProcessor = require("./molecular-processor");
const { resolveName, getPropertiesByCID } = require("./name-resolver");

// Text-first API surface (UI prioritizes text input)
class Structuralizer {
  constructor(apiKey) {
    this.isTestMode = process.env.NODE_ENV === 'test';
    this.client = OpenAIClient ? new OpenAIClient({ apiKey: apiKey || process.env.OPENAI_API_KEY || '' }) : null;
    this.isOpenAIAvailable = !!this.client && !!(apiKey || process.env.OPENAI_API_KEY);
    this.modelName = process.env.OPENAI_MODEL || process.env.OPENAI_DEFAULT_MODEL || 'gpt-4o';
    this.chemicalInstructions = imageToSmiles_instructions();
    this.molecularProcessor = new MolecularProcessor();
  }

  // Unified multimodal entry: accepts { object?, imageBase64? } and returns compounds
  async analyze(payload) {
    const object = typeof payload?.object === 'string' ? payload.object : '';
    const imageBase64 = typeof payload?.imageBase64 === 'string' ? payload.imageBase64 : null;
    let mergedNames = [];
    let mergedChemicals = [];

    // Text path (preferred)
    if (object && object.trim()) {
      try {
        const byText = await this.analyzeText(object.trim());
        const list = Array.isArray(byText?.molecules) ? byText.molecules : [];
        mergedNames.push(...list.map(m => ({ name: m.name })));
        const chems = Array.isArray(byText?.chemicals) ? byText.chemicals : [];
        mergedChemicals.push(...chems);
      } catch (_) {}
    }

    // Optional image path
    if (imageBase64 && this.isOpenAIAvailable) {
      try {
        const byImage = await this.analyzeImage(imageBase64);
        const chems = Array.isArray(byImage?.chemicals) ? byImage.chemicals : [];
        mergedChemicals.push(...chems);
        mergedNames.push(...chems.filter(c => c?.name).map(c => ({ name: c.name })));
      } catch (_) {}
    }

    // Deduplicate any immediate chemicals
    const byKey = new Map();
    for (const c of mergedChemicals) {
      const key = (c.sdfPath || c.smiles || c.name || '').toLowerCase();
      if (key && !byKey.has(key)) byKey.set(key, c);
    }
    const chemicals = Array.from(byKey.values());
    if (chemicals.length > 0) {
      return { object: object || (imageBase64 ? 'image' : ''), chemicals };
    }

    // Resolve names → structures
    const nameSet = new Set(mergedNames.map(n => (n?.name || '').toLowerCase()).filter(Boolean));
    const uniqueNames = Array.from(nameSet);
    const resolved = [];
    for (const nm of uniqueNames) {
      try {
        const byName = await this.molecularProcessor.generateSDFByName(nm, false);
        if (byName && byName.sdfPath) resolved.push({ name: byName.name || nm, sdfPath: byName.sdfPath });
      } catch (_) {}
    }
    return { object: object || (imageBase64 ? 'image' : ''), chemicals: resolved };
  }

  // Primary: text → names → structures (UI default)
  async analyzeText(object) {
    if (!this.isOpenAIAvailable) throw new Error("AI service unavailable for names-only extraction");
    if (this.isTestMode) {
      const namesPayload = await this.extractNamesOnly(object);
      const list = Array.isArray(namesPayload?.molecules) ? namesPayload.molecules : [];
      const smilesPayload = await this.convertNamesToSmilesProgrammatically({ object, molecules: list });
      const chemicals = Array.isArray(smilesPayload?.molecules) ? smilesPayload.molecules : [];
      return { object: namesPayload.object || object, chemicals, molecules: list };
    }
    const namesPayload = await this.extractNamesOnly(object);
    const namesList = Array.isArray(namesPayload?.molecules) ? namesPayload.molecules : [];
    if (namesList.length === 0) return { object, chemicals: [], molecules: [], meta: { strategy: 'two-step', names: [], namesCount: 0, sdfCount: 0 } };
    const resolved = [];
    for (const m of namesList) {
      const res = await this.molecularProcessor.generateSDFByName(m.name, false);
      if (res && res.sdfPath) resolved.push({ name: res.name || m.name, sdfPath: res.sdfPath });
    }
    const names = namesList.map(n => n.name);
    const meta = { strategy: 'two-step', names, namesCount: names.length, sdfCount: resolved.length };
    return { object, chemicals: resolved, molecules: resolved, meta };
  }

  // Secondary: image analysis (used from camera UI)
  async analyzeImage(imageBase64, croppedImageBase64 = null, x = null, y = null, cropMiddleX = null, cropMiddleY = null, cropSize = null) {
    if (!imageBase64 || imageBase64.trim().length === 0) throw new Error("Empty image data");
    if (!this.isOpenAIAvailable) throw new Error("AI service unavailable");
    if (/^https?:\/\//i.test(imageBase64)) {
      const fetch = require('node-fetch');
      const resp = await fetch(imageBase64);
      if (!resp.ok) throw new Error(`Failed to fetch image URL: ${resp.status}`);
      const buf = await resp.buffer();
      imageBase64 = buf.toString('base64');
    }
    const messages = [{ role: "user", content: [{ type: "text", text: this.chemicalInstructions }, { type: "image_url", image_url: { url: `data:image/jpeg;base64,${imageBase64}`, detail: "high" } }] }];
    if (croppedImageBase64) {
      messages[0].content.push({ type: "text", text: `Here's a cropped view of the area of interest. Analyze the chemical composition of the material or substance visible in this region. Use the examples above as your guide for accurate SMILES notation:` });
      messages[0].content.push({ type: "image_url", image_url: { url: `data:image/jpeg;base64,${croppedImageBase64}`, detail: "high" } });
    }
    const response = await this.client.chat.completions.create({ model: this.modelName, messages, max_tokens: 1000, temperature: 0.1, response_format: { type: "json_object" } });
    const content = response.choices[0].message.content;
    const parsed = JSON.parse(content);
    return { object: parsed.object || "Unknown object", chemicals: parsed.chemicals || [] };
  }

  async extractNamesOnly(object) {
    if (!this.isOpenAIAvailable) return { object, molecules: [] };
    const prompt = textToNames_prompt(object);
    const response = await this.client.chat.completions.create({ model: this.modelName, messages: [{ role: "user", content: prompt }], max_tokens: 800, temperature: 0.1, response_format: { type: "json_object" } });
    const content = response.choices[0].message.content;
    let parsed;
    try { parsed = JSON.parse(content); } catch (_) { const match = content && content.match(/\{[\s\S]*\}/); parsed = match ? JSON.parse(match[0]) : { object, molecules: [] }; }
    let molecules = [];
    if (Array.isArray(parsed.molecules)) molecules = parsed.molecules; else if (Array.isArray(parsed.chemicals)) molecules = parsed.chemicals.map((c) => ({ name: c?.name || c?.title || c?.iupac || '', cid: c?.cid ?? null })).filter((m) => typeof m.name === 'string' && m.name.trim().length > 0);
    return { object: parsed.object || object, molecules };
  }

  async convertNamesToSmilesProgrammatically(namesPayload) {
    const object = namesPayload?.object || "";
    const input = Array.isArray(namesPayload?.molecules) ? namesPayload.molecules : [];
    const results = [];
    for (const mol of input) {
      const name = mol?.name || '';
      let cid = mol?.cid ?? null;
      let smiles = null;
      try {
        if (cid) {
          const props = await getPropertiesByCID(cid);
          smiles = props?.smiles || null;
          if (!smiles) { const res = await resolveName(name); cid = res?.cid || cid; smiles = res?.smiles || null; }
        } else {
          const res = await resolveName(name); cid = res?.cid || null; smiles = res?.smiles || null;
        }
      } catch (_) {}
      results.push({ name, cid, smiles: smiles || null, status: smiles ? 'ok' : 'lookup_required' });
    }
    return { object, molecules: results };
  }

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


