/**
 * CACHE WRAPPER FOR ASYNC FUNCTIONS
 * ==================================
 * User prompt (from Structuralizer.js lines 102-119, 291-317):
 * "extract the cache components into it's own file, using the simplest way to integrate
 * with existing code when needed"
 * 
 * DESIGN RATIONALE:
 * - Wraps expensive async operations with caching
 * - Extracted hash function from Structuralizer._hashString (lines 309-317)
 * - Extracted cache key generation from Structuralizer._getCacheKey (lines 291-303)
 * - Can work with any cache implementing get(key) / set(key, value, ttl)
 * 
 * DECORATOR CONSIDERATION:
 * User asked about using decorators/attributes. JavaScript decorators are still Stage 3 proposal.
 * Using them would require:
 * - @babel/plugin-proposal-decorators (adds build complexity)
 * - OR experimental Node.js flags (unstable)
 * - Less clear for developers unfamiliar with decorators
 * 
 * CONCLUSION: Simple function wrapper is better - no build tools, works everywhere, clearer intent
 * 
 * USAGE EXAMPLES:
 * 
 * @example Basic usage
 * const wrapper = new CacheWrapper({ cache: myCache, enabled: true });
 * 
 * const result = await wrapper.wrap(
 *   'my-operation',
 *   { userId: 123, query: 'data' },
 *   async () => expensiveOperation(),
 *   { ttl: 60000 }
 * );
 * 
 * @example With hash for large payloads
 * const result = await wrapper.wrap(
 *   'image-analysis',
 *   {
 *     object: 'coffee',
 *     imageHash: CacheWrapper.hashString(largeImageBase64.substring(0, 100))
 *   },
 *   async () => analyzeImage(largeImageBase64),
 *   { ttl: 300000 }
 * );
 */

class CacheWrapper {
  constructor(options = {}) {
    this.cache = options.cache || null;
    this.enabled = options.enabled !== false;
    this.logger = options.logger || console;
    this.defaultTTL = options.defaultTTL || 300000; // 5 minutes
  }

  /**
   * Wrap an async function with caching
   * Integrates cache logic from Structuralizer.chemicals (lines 102-119)
   * 
   * @param {string} prefix - Cache key prefix (e.g., 'structuralizer')
   * @param {Object} payload - Data to include in cache key
   * @param {Function} fn - Async function to execute on cache miss
   * @param {Object} options - { ttl: number }
   * @returns {Promise<any>}
   */
  async wrap(prefix, payload, fn, options = {}) {
    const startTime = Date.now();
    const ttl = options.ttl || this.defaultTTL;
    
    // Bypass cache if disabled or not available
    if (!this.cache || !this.enabled) {
      return await fn();
    }

    // Generate cache key
    const cacheKey = this._generateKey(prefix, payload);
    
    // Try to get from cache
    const cached = await this.cache.get(cacheKey);
    if (cached) {
      this.logger.info?.('Cache hit', {
        prefix,
        duration: Date.now() - startTime
      });
      return cached;
    }
    
    // Cache miss - execute function
    const result = await fn();
    
    // Store in cache if result is valid
    if (result !== null && result !== undefined) {
      await this.cache.set(cacheKey, result, ttl);
    }
    
    this.logger.info?.('Cache miss', {
      prefix,
      duration: Date.now() - startTime
    });
    
    return result;
  }

  /**
   * Wrap with a pre-generated cache key
   * @param {string} cacheKey - Full cache key
   * @param {Function} fn - Async function to execute
   * @param {Object} options - { ttl: number }
   * @returns {Promise<any>}
   */
  async wrapWithKey(cacheKey, fn, options = {}) {
    const ttl = options.ttl || this.defaultTTL;
    
    if (!this.cache || !this.enabled) {
      return await fn();
    }

    const cached = await this.cache.get(cacheKey);
    if (cached) {
      return cached;
    }
    
    const result = await fn();
    
    if (result !== null && result !== undefined) {
      await this.cache.set(cacheKey, result, ttl);
    }
    
    return result;
  }

  /**
   * Generate cache key from prefix and payload
   * Extracted from Structuralizer._getCacheKey (lines 291-303)
   * 
   * @param {string} prefix
   * @param {Object} payload
   * @returns {string}
   */
  _generateKey(prefix, payload) {
    return `${prefix}:${JSON.stringify(payload)}`;
  }

  /**
   * Simple string hash for cache keys
   * Extracted from Structuralizer._hashString (lines 309-317)
   * 
   * USAGE: Use this to hash large strings (like images) before including in cache key
   * 
   * @param {string} str - String to hash
   * @returns {string} - Base36 hash
   * 
   * @example
   * const imageHash = CacheWrapper.hashString(imageBase64.substring(0, 100));
   * await wrapper.wrap('detection', { object: 'coffee', imageHash }, () => detect());
   */
  static hashString(str) {
    let hash = 0;
    for (let i = 0; i < str.length; i++) {
      const char = str.charCodeAt(i);
      hash = ((hash << 5) - hash) + char;
      hash = hash & hash; // Convert to 32-bit integer
    }
    return hash.toString(36);
  }
}

module.exports = CacheWrapper;
