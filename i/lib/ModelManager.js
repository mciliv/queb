const { report } = require('./apm/client');

class ModelManager {
  constructor(aiService) {
    this.aiService = aiService;
    this.status = new Map();
    this.backoffTime = 60000;
  }

  async callWithFailover(params, priorityList) {
    let lastError = null;

    for (const provider of priorityList) {
      if (this._isThrottled(provider)) continue;

      try {
        this.aiService.switchProvider(provider);
        return await this.aiService.callAPI(params);
      } catch (error) {
        lastError = error;

        if (this._isRateLimitError(error)) {
          report({ status: 'failover', provider, error: '429 rate limit' });
          this._markThrottled(provider);
          continue;
        }

        console.error(`[ModelManager] Provider ${provider} failed:`, error.message);
        continue;
      }
    }

    throw new Error(`All providers failed. Last error: ${lastError?.message}`);
  }

  _isRateLimitError(error) {
    const msg = error.message.toLowerCase();
    return msg.includes('429') || msg.includes('rate limit') || msg.includes('too many requests');
  }

  _markThrottled(provider) {
    this.status.set(provider, {
      throttledUntil: Date.now() + this.backoffTime
    });
  }

  _isThrottled(provider) {
    const state = this.status.get(provider);
    if (!state) return false;
    return Date.now() < state.throttledUntil;
  }
}

module.exports = ModelManager;
