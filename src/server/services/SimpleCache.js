/**
 * SIMPLE IN-MEMORY CACHE WITH TTL
 * ================================
 * User prompt (from Structuralizer.js lines 102-119, 291-317):
 * "extract the cache components into it's own file, using the simplest way to integrate 
 * with existing code when needed, but for now simply describe how to use it"
 * 
 * DESIGN RATIONALE:
 * - No external dependencies (no Redis, no node-cache)
 * - Simple Map-based implementation with TTL
 * - Compatible with any cache interface expecting get(key) / set(key, value, ttl)
 * - Can be replaced with Redis/Memcached later without changing consumer code
 * 
 * STANDALONE POTENTIAL:
 * This can be extracted as @queb/simple-cache or simple-ttl-cache
 * It's small (~80 lines), has no dependencies, and solves a common problem
 */

class SimpleCache {
  constructor(options = {}) {
    this.maxSize = options.maxSize || 1000;
    this.defaultTTL = options.defaultTTL || 300000; // 5 minutes
    this.cleanupInterval = options.cleanupInterval || 60000; // 1 minute
    
    // Storage: Map<key, { value, expiresAt }>
    this.store = new Map();
    
    // Stats for monitoring
    this.stats = {
      hits: 0,
      misses: 0,
      sets: 0,
      evictions: 0
    };
    
    // Start cleanup timer
    this._startCleanup();
  }
  
  /**
   * Get value from cache
   * @param {string} key
   * @returns {Promise<any|null>}
   */
  async get(key) {
    const entry = this.store.get(key);
    
    if (!entry) {
      this.stats.misses++;
      return null;
    }
    
    // Check if expired
    if (Date.now() > entry.expiresAt) {
      this.store.delete(key);
      this.stats.misses++;
      return null;
    }
    
    this.stats.hits++;
    return entry.value;
  }
  
  /**
   * Set value in cache
   * @param {string} key
   * @param {any} value
   * @param {number} ttl - Time to live in milliseconds
   * @returns {Promise<void>}
   */
  async set(key, value, ttl) {
    const timeToLive = ttl || this.defaultTTL;
    
    // Evict oldest if at capacity
    if (this.store.size >= this.maxSize && !this.store.has(key)) {
      const firstKey = this.store.keys().next().value;
      this.store.delete(firstKey);
      this.stats.evictions++;
    }
    
    this.store.set(key, {
      value,
      expiresAt: Date.now() + timeToLive
    });
    
    this.stats.sets++;
  }
  
  /**
   * Delete key from cache
   * @param {string} key
   * @returns {Promise<boolean>}
   */
  async delete(key) {
    return this.store.delete(key);
  }
  
  /**
   * Clear all cache
   * @returns {Promise<void>}
   */
  async clear() {
    this.store.clear();
    this.stats = {
      hits: 0,
      misses: 0,
      sets: 0,
      evictions: 0
    };
  }
  
  /**
   * Get cache statistics
   * @returns {Object}
   */
  getStats() {
    return {
      ...this.stats,
      size: this.store.size,
      hitRate: this.stats.hits + this.stats.misses > 0
        ? (this.stats.hits / (this.stats.hits + this.stats.misses) * 100).toFixed(2) + '%'
        : '0%'
    };
  }
  
  /**
   * Start cleanup timer to remove expired entries
   * @private
   */
  _startCleanup() {
    this.cleanupTimer = setInterval(() => {
      const now = Date.now();
      let cleaned = 0;
      
      for (const [key, entry] of this.store.entries()) {
        if (now > entry.expiresAt) {
          this.store.delete(key);
          cleaned++;
        }
      }
      
      if (cleaned > 0 && process.env.NODE_ENV === 'development') {
        console.log(`[SimpleCache] Cleaned up ${cleaned} expired entries`);
      }
    }, this.cleanupInterval);
    
    // Don't keep process alive
    this.cleanupTimer.unref();
  }
  
  /**
   * Stop cleanup timer (for graceful shutdown)
   */
  destroy() {
    if (this.cleanupTimer) {
      clearInterval(this.cleanupTimer);
    }
  }
}

module.exports = SimpleCache;
