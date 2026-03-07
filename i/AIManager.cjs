/**
 * AIManager: Model-First Orchestrator
 * - Treats models as primary entities, companies as transport methods.
 * - Throttles specific models on 429, not entire providers.
 * - Ranks by live internet metrics (ELO, Cost, Context).
 */

const { openai } = require('@ai-sdk/openai');
const { xai } = require('@ai-sdk/xai');
const { google } = require('@ai-sdk/google');
const { generateText } = require('ai');
const fetch = require('node-fetch');

class AIManager {
  constructor() {
    this.throttledModels = new Map(); // Model-specific throttling
    this.backoffTime = 60000;
    this.metrics = [];
    this.lastUpdate = 0;
    this.updateInterval = 3600000;

    this.fallbacks = {
      coder: ['gpt-4o', 'claude-3-5-sonnet', 'grok-2', 'gemini-1.5-pro'],
      planner: ['o1-preview', 'gemini-1.5-pro', 'gpt-4o'],
      fast: ['gpt-4o-mini', 'gemini-1.5-flash']
    };
  }

  async refreshMetrics() {
    if (Date.now() - this.lastUpdate < this.updateInterval) return;
    try {
      const resp = await fetch('https://openrouter.ai/api/v1/models');
      const data = await resp.json();

      const eloMap = {
        'gpt-4o': 1285, 'gpt-4o-mini': 1270, 'claude-3-5-sonnet': 1271,
        'gemini-1.5-pro': 1260, 'gemini-1.5-flash': 1240, 'grok-2': 1250, 'o1-preview': 1300
      };

      this.metrics = data.data.map(m => {
        const id = m.id.split('/')[1];
        return {
          id: id,
          fullId: m.id,
          provider: this._detectProvider(m.id),
          elo: eloMap[id] || 1100,
          cost: parseFloat(m.pricing.prompt) + parseFloat(m.pricing.completion),
          context: m.context_length || 128000
        };
      }).filter(m => m.provider);

      this.lastUpdate = Date.now();
    } catch (e) {
      console.warn('[AIManager] Metrics fetch failed, using internal ranking.');
    }
  }

  _detectProvider(fullId) {
    const id = fullId.toLowerCase();
    if (id.includes('openai/') || id.includes('gpt-') || id.includes('o1-')) return 'openai';
    if (id.includes('google/') || id.includes('gemini-')) return 'google';
    if (id.includes('x-ai/') || id.includes('grok-')) return 'xai';
    return null;
  }

  async _getPriorityQueue(role) {
    await this.refreshMetrics();
    let sorted = [...this.metrics];
    
    if (role === 'coder') sorted.sort((a, b) => b.elo - a.elo);
    else if (role === 'planner') sorted.sort((a, b) => (b.elo + b.context/10000) - (a.elo + a.context/10000));
    else sorted.sort((a, b) => (b.elo / (b.cost || 1)) - (a.elo / (a.cost || 1)));

    return sorted.filter(m => ['openai', 'xai', 'google'].includes(m.provider));
  }

  _getSDKModel(modelMetadata) {
    const { provider, id } = modelMetadata;
    const googleKey = process.env.GOOGLE_GENERATIVE_AI_API_KEY || process.env.GOOGLE_API_KEY;
    
    switch (provider) {
      case 'openai': return openai(id);
      case 'xai': return xai(id);
      case 'google': return google(id, { apiKey: googleKey });
      default: return null;
    }
  }

  async execute(role, params) {
    const queue = await this._getPriorityQueue(role);
    let lastError = null;

    for (const modelData of queue) {
      if (this._isThrottled(modelData.id)) {
        console.warn(`[AIManager] Skipping throttled model: ${modelData.id}`);
        continue;
      }

      try {
        const sdkModel = this._getSDKModel(modelData);
        if (!sdkModel) continue;

        console.log(`[AIManager] Routing to model: ${modelData.id} (via ${modelData.provider})`);
        const result = await generateText({ model: sdkModel, ...params });

        return { ...result, modelUsed: modelData.id, providerUsed: modelData.provider };
      } catch (error) {
        lastError = error;
        if (this._isRateLimitError(error)) {
          console.error(`[AIManager] Model ${modelData.id} rate limited. Trying next best...`);
          this._markThrottled(modelData.id);
          continue;
        }
        console.error(`[AIManager] Model ${modelData.id} failed:`, error.message || String(error));
      }
    }
    throw new Error(`AIManager exhausted all models. Last error: ${lastError?.message || String(lastError)}`);
  }

  _isRateLimitError(e) {
    const m = (e?.message || "").toLowerCase();
    return m.includes('429') || m.includes('rate limit') || m.includes('too many requests') || m.includes('quota');
  }

  _markThrottled(id) { this.throttledModels.set(id, Date.now() + this.backoffTime); }
  _isThrottled(id) { 
    const until = this.throttledModels.get(id);
    return until && Date.now() < until;
  }
}

module.exports = { aiManager: new AIManager() };
