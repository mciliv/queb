const { fetchHighestQualityIngredients } = require('./ingredient-delivery');

class DeliveryService {
  constructor({ aiService, promptEngine, logger }) {
    this.aiService = aiService;
    this.promptEngine = promptEngine;
    this.logger = logger;
  }

  async analyzeDeliveryOptions({ item, city, lat, lng }) {
    if (!item || typeof item !== 'string') {
      throw new Error('Item name is required');
    }

    const [aiAnalysis, supplierResults] = await Promise.allSettled([
      this._getAIAnalysis(item, { city, lat, lng }),
      this._getSupplierPricing({ city, lat, lng })
    ]);

    const analysis = aiAnalysis.status === 'fulfilled' ? aiAnalysis.value : null;
    const suppliers = supplierResults.status === 'fulfilled' ? supplierResults.value : [];

    if (aiAnalysis.status === 'rejected') {
      this.logger.warn('AI delivery analysis failed', { error: aiAnalysis.reason.message });
    }
    if (supplierResults.status === 'rejected') {
      this.logger.warn('Supplier lookup failed', { error: supplierResults.reason.message });
    }

    return this._mergeResults(item, analysis, suppliers, { city, lat, lng });
  }

  async _getAIAnalysis(item, location) {
    const prompt = this.promptEngine.generateDeliveryPrompt(item, location);

    const result = await this.aiService.callAPI({
      messages: [
        { role: 'system', content: 'You are a procurement and sourcing analyst API. Always respond with valid JSON only, no markdown, no commentary.' },
        { role: 'user', content: prompt }
      ],
      max_tokens: 4000,
      temperature: 1.0
    });

    if (!result || !result.item) {
      throw new Error('Invalid delivery analysis response from AI');
    }

    return result;
  }

  async _getSupplierPricing({ city, lat, lng }) {
    if (!city && (lat === undefined || lng === undefined)) {
      return [];
    }

    try {
      const results = await fetchHighestQualityIngredients({
        city, lat, lng, limit: 10, timeoutMs: 8000
      });
      return Array.isArray(results) ? results : [];
    } catch (err) {
      this.logger.warn('Supplier pricing unavailable', { error: err.message });
      return [];
    }
  }

  _mergeResults(item, aiAnalysis, suppliers, location) {
    const base = aiAnalysis || {
      item,
      growable: { possible: false, difficulty: 'unknown', notes: 'AI analysis unavailable' },
      purchase_options: [],
      recommended_strategy: 'buy',
      reasoning: 'Defaulting to purchase - AI analysis was unavailable',
      time_to_obtain: { buy_days: 2, grow_days: null }
    };

    // Merge local supplier data into purchase options
    const supplierOptions = suppliers.map(s => ({
      source_type: 'local_store',
      vendor_name: s.supplier || s.name,
      estimated_price_usd: null,
      unit: 'each',
      availability: s.freshness_hours <= 12 ? 'immediate' : 'seasonal',
      quality_tier: s.quality_score >= 95 ? 'premium' : s.quality_score >= 90 ? 'standard' : 'bulk',
      quality_score: s.quality_score,
      supplier_rating: s.supplier_rating,
      distance_km: s.distance_km,
      item_name: s.name
    }));

    const allOptions = [
      ...(base.purchase_options || []),
      ...supplierOptions
    ];

    return {
      item: base.item || item,
      location: {
        city: location.city || null,
        lat: location.lat || null,
        lng: location.lng || null
      },
      growable: base.growable,
      purchase_options: allOptions,
      local_suppliers: supplierOptions,
      recommended_strategy: base.recommended_strategy,
      reasoning: base.reasoning,
      time_to_obtain: base.time_to_obtain
    };
  }
}

module.exports = DeliveryService;
