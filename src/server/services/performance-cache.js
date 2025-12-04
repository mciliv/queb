/**
 * Performance Cache Service
 * 
 * Combines performance tracking and intelligent caching using MCP integration.
 * Simple to use, provides immediate value.
 * 
 * Usage:
 *   const cache = new PerformanceCache();
 *   const result = await cache.wrap('analyzeText', input, async () => {
 *     return await expensiveOperation(input);
 *   });
 */

const crypto = require('crypto');
const { getMCPService } = require('./mcp-integration');

class PerformanceCache {
  constructor(options = {}) {
    this.mcpService = getMCPService();
    this.logger = options.logger || console;
    this.defaultTTL = options.ttl || 7 * 24 * 60 * 60 * 1000; // 7 days
    this.enablePerformanceTracking = options.trackPerformance !== false;
    this.enableCaching = options.enableCache !== false;
  }

  /**
   * Wrap an expensive operation with caching and performance tracking
   * 
   * @param {string} operation - Name of the operation (e.g., 'analyzeText')
   * @param {any} input - Input that determines cache key
   * @param {Function} fn - Async function to execute if cache miss
   * @param {Object} options - { ttl, forceRefresh, metadata }
   * @returns {Promise<any>} - Result from cache or function execution
   */
  async wrap(operation, input, fn, options = {}) {
    const startTime = Date.now();
    const cacheKey = this._generateCacheKey(operation, input);
    
    try {
      // Check cache first (unless force refresh)
      if (this.enableCaching && !options.forceRefresh) {
        const cached = await this._getFromCache(operation, cacheKey);
        if (cached !== null) {
          const duration = Date.now() - startTime;
          
          if (this.enablePerformanceTracking) {
            await this._trackPerformance(operation, duration, {
              ...options.metadata,
              cacheHit: true,
              inputSize: this._estimateSize(input)
            });
          }
          
          this.logger.info(`✅ Cache hit for ${operation} (${duration}ms)`);
          return cached;
        }
      }
      
      // Execute function
      this.logger.info(`🔄 Cache miss for ${operation}, executing...`);
      const result = await fn();
      
      const duration = Date.now() - startTime;
      
      // Store in cache
      if (this.enableCaching && result !== null && result !== undefined) {
        await this._setInCache(operation, cacheKey, input, result, options.ttl);
      }
      
      // Track performance
      if (this.enablePerformanceTracking) {
        await this._trackPerformance(operation, duration, {
          ...options.metadata,
          cacheHit: false,
          inputSize: this._estimateSize(input),
          resultSize: this._estimateSize(result)
        });
      }
      
      this.logger.info(`✅ ${operation} completed (${duration}ms)`);
      return result;
      
    } catch (error) {
      const duration = Date.now() - startTime;
      
      // Track failed operation
      if (this.enablePerformanceTracking) {
        await this._trackPerformance(operation, duration, {
          ...options.metadata,
          success: false,
          error: error.message
        });
      }
      
      throw error;
    }
  }

  /**
   * Get performance statistics for an operation
   */
  async getStats(operation, hours = 24) {
    if (!this.enablePerformanceTracking) {
      return null;
    }

    try {
      const stats = await this.mcpService.executeReadQuery(`
        SELECT 
          operation,
          COUNT(*) as total_calls,
          COUNT(*) FILTER (WHERE (metadata->>'cacheHit')::boolean = true) as cache_hits,
          COUNT(*) FILTER (WHERE (metadata->>'cacheHit')::boolean = false) as cache_misses,
          AVG(duration_ms) as avg_duration_ms,
          MIN(duration_ms) as min_duration_ms,
          MAX(duration_ms) as max_duration_ms,
          AVG(CASE WHEN (metadata->>'cacheHit')::boolean = false THEN duration_ms END) as avg_uncached_ms
        FROM performance_metrics
        WHERE operation = $1 
          AND timestamp > NOW() - INTERVAL '${hours} hours'
        GROUP BY operation
      `, [operation]);

      if (stats.length === 0) {
        return null;
      }

      const stat = stats[0];
      return {
        operation: stat.operation,
        totalCalls: parseInt(stat.total_calls),
        cacheHits: parseInt(stat.cache_hits || 0),
        cacheMisses: parseInt(stat.cache_misses || 0),
        cacheHitRate: stat.cache_hits 
          ? (parseInt(stat.cache_hits) / parseInt(stat.total_calls) * 100).toFixed(1) + '%'
          : '0%',
        avgDuration: Math.round(stat.avg_duration_ms) + 'ms',
        minDuration: Math.round(stat.min_duration_ms) + 'ms',
        maxDuration: Math.round(stat.max_duration_ms) + 'ms',
        avgUncached: stat.avg_uncached_ms ? Math.round(stat.avg_uncached_ms) + 'ms' : 'N/A'
      };
    } catch (error) {
      this.logger.warn('Could not get stats (table may not exist):', error.message);
      return null;
    }
  }

  /**
   * Clear cache for a specific operation or all
   */
  async clearCache(operation = null) {
    try {
      if (operation) {
        await this.mcpService.executeReadQuery(
          'DELETE FROM analysis_cache WHERE operation = $1',
          [operation]
        );
        this.logger.info(`🗑️  Cleared cache for ${operation}`);
      } else {
        await this.mcpService.executeReadQuery('DELETE FROM analysis_cache');
        this.logger.info('🗑️  Cleared all cache');
      }
    } catch (error) {
      this.logger.warn('Could not clear cache:', error.message);
    }
  }

  /**
   * Get overall cache statistics
   */
  async getCacheStats() {
    try {
      const stats = await this.mcpService.executeReadQuery(`
        SELECT 
          operation,
          COUNT(*) as entries,
          SUM(hit_count) as total_hits,
          AVG(hit_count) as avg_hits,
          MAX(created_at) as last_used
        FROM analysis_cache
        GROUP BY operation
        ORDER BY total_hits DESC
      `);

      return stats.map(s => ({
        operation: s.operation,
        entries: parseInt(s.entries),
        totalHits: parseInt(s.total_hits),
        avgHits: parseFloat(s.avg_hits).toFixed(1),
        lastUsed: s.last_used
      }));
    } catch (error) {
      this.logger.warn('Could not get cache stats:', error.message);
      return [];
    }
  }

  /**
   * PRIVATE METHODS
   */

  async _getFromCache(operation, cacheKey) {
    try {
      const cached = await this.mcpService.executeReadQuery(`
        SELECT result, created_at, hit_count
        FROM analysis_cache
        WHERE cache_key = $1 
          AND operation = $2
          AND created_at > NOW() - INTERVAL '7 days'
        LIMIT 1
      `, [cacheKey, operation]);

      if (cached.length > 0) {
        // Increment hit count asynchronously
        this._incrementHitCount(cacheKey, operation).catch(() => {});
        
        return JSON.parse(cached[0].result);
      }
    } catch (error) {
      // Table might not exist, silent fail
      this.logger.debug('Cache read failed:', error.message);
    }
    
    return null;
  }

  async _setInCache(operation, cacheKey, input, result, ttl = null) {
    try {
      const inputText = typeof input === 'string' ? input : JSON.stringify(input);
      
      await this.mcpService.executeReadQuery(`
        INSERT INTO analysis_cache 
        (operation, cache_key, input_text, result, created_at, hit_count)
        VALUES ($1, $2, $3, $4, NOW(), 0)
        ON CONFLICT (cache_key, operation) 
        DO UPDATE SET 
          result = $4, 
          created_at = NOW(),
          hit_count = 0
      `, [
        operation,
        cacheKey,
        inputText.substring(0, 1000), // Limit input text storage
        JSON.stringify(result)
      ]);
    } catch (error) {
      // Silent fail - caching is not critical
      this.logger.debug('Cache write failed:', error.message);
    }
  }

  async _incrementHitCount(cacheKey, operation) {
    try {
      await this.mcpService.executeReadQuery(`
        UPDATE analysis_cache 
        SET hit_count = hit_count + 1 
        WHERE cache_key = $1 AND operation = $2
      `, [cacheKey, operation]);
    } catch (error) {
      // Silent fail
    }
  }

  async _trackPerformance(operation, duration, metadata = {}) {
    try {
      await this.mcpService.executeReadQuery(`
        INSERT INTO performance_metrics 
        (operation, duration_ms, metadata, timestamp)
        VALUES ($1, $2, $3, NOW())
      `, [operation, duration, JSON.stringify(metadata)]);
    } catch (error) {
      // Silent fail - tracking is not critical
      this.logger.debug('Performance tracking failed:', error.message);
    }
  }

  _generateCacheKey(operation, input) {
    const normalized = typeof input === 'string' 
      ? input.toLowerCase().trim() 
      : JSON.stringify(input);
    
    return crypto
      .createHash('sha256')
      .update(operation + ':' + normalized)
      .digest('hex');
  }

  _estimateSize(obj) {
    if (typeof obj === 'string') {
      return obj.length;
    }
    return JSON.stringify(obj).length;
  }
}

module.exports = PerformanceCache;









