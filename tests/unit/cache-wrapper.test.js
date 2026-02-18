/**
 * CACHE WRAPPER TESTS
 * ===================
 * User prompt: "extract the cache components into it's own file"
 * 
 * Tests for SimpleCache and CacheWrapper extracted from Structuralizer
 */

const SimpleCache = require('../../src/server/services/SimpleCache');
const CacheWrapper = require('../../src/server/services/CacheWrapper');

describe('SimpleCache', () => {
  let cache;

  beforeEach(() => {
    cache = new SimpleCache({
      maxSize: 10,
      defaultTTL: 1000,
      cleanupInterval: 5000
    });
  });

  afterEach(() => {
    cache.destroy();
  });

  test('stores and retrieves values', async () => {
    await cache.set('key1', 'value1', 1000);
    const value = await cache.get('key1');
    expect(value).toBe('value1');
  });

  test('returns null for non-existent keys', async () => {
    const value = await cache.get('nonexistent');
    expect(value).toBeNull();
  });

  test('respects TTL and expires entries', async () => {
    await cache.set('key1', 'value1', 100);
    
    // Immediate read - should exist
    const value1 = await cache.get('key1');
    expect(value1).toBe('value1');
    
    // Wait for expiration
    await new Promise(resolve => setTimeout(resolve, 150));
    
    // Should be expired
    const value2 = await cache.get('key1');
    expect(value2).toBeNull();
  });

  test('tracks cache statistics', async () => {
    await cache.set('key1', 'value1', 1000);
    await cache.get('key1'); // hit
    await cache.get('key2'); // miss
    
    const stats = cache.getStats();
    expect(stats.hits).toBe(1);
    expect(stats.misses).toBe(1);
    expect(stats.sets).toBe(1);
    expect(stats.size).toBe(1);
  });

  test('evicts oldest entry when at capacity', async () => {
    // Fill to capacity
    for (let i = 0; i < 10; i++) {
      await cache.set(`key${i}`, `value${i}`, 10000);
    }
    
    expect(cache.getStats().size).toBe(10);
    
    // Add one more - should evict oldest
    await cache.set('key10', 'value10', 10000);
    
    const stats = cache.getStats();
    expect(stats.size).toBe(10);
    expect(stats.evictions).toBe(1);
  });

  test('clears all entries', async () => {
    await cache.set('key1', 'value1', 1000);
    await cache.set('key2', 'value2', 1000);
    
    await cache.clear();
    
    const stats = cache.getStats();
    expect(stats.size).toBe(0);
    expect(await cache.get('key1')).toBeNull();
  });
});

describe('CacheWrapper', () => {
  let cache;
  let wrapper;
  let logger;

  beforeEach(() => {
    cache = new SimpleCache({ maxSize: 100 });
    logger = {
      info: jest.fn()
    };
    wrapper = new CacheWrapper({
      cache,
      enabled: true,
      logger
    });
  });

  afterEach(() => {
    cache.destroy();
  });

  test('caches expensive operations', async () => {
    let callCount = 0;
    const expensiveOp = async () => {
      callCount++;
      return { result: 'expensive data' };
    };

    // First call - cache miss
    const result1 = await wrapper.wrap(
      'operation',
      { id: 1 },
      expensiveOp,
      { ttl: 1000 }
    );
    expect(callCount).toBe(1);
    expect(result1).toEqual({ result: 'expensive data' });

    // Second call - cache hit
    const result2 = await wrapper.wrap(
      'operation',
      { id: 1 },
      expensiveOp,
      { ttl: 1000 }
    );
    expect(callCount).toBe(1); // Not called again
    expect(result2).toEqual(result1);
  });

  test('different payloads create different cache keys', async () => {
    let callCount = 0;
    const expensiveOp = async (n) => {
      callCount++;
      return { result: n };
    };

    await wrapper.wrap('op', { id: 1 }, () => expensiveOp(1));
    await wrapper.wrap('op', { id: 2 }, () => expensiveOp(2));

    // Both should be called - different cache keys
    expect(callCount).toBe(2);
  });

  test('bypasses cache when disabled', async () => {
    const disabledWrapper = new CacheWrapper({
      cache,
      enabled: false,
      logger
    });

    let callCount = 0;
    const expensiveOp = async () => {
      callCount++;
      return { result: 'data' };
    };

    await disabledWrapper.wrap('op', { id: 1 }, expensiveOp);
    await disabledWrapper.wrap('op', { id: 1 }, expensiveOp);

    // Both calls should execute
    expect(callCount).toBe(2);
  });

  test('bypasses cache when cache is null', async () => {
    const nullCacheWrapper = new CacheWrapper({
      cache: null,
      enabled: true,
      logger
    });

    let callCount = 0;
    const expensiveOp = async () => {
      callCount++;
      return { result: 'data' };
    };

    await nullCacheWrapper.wrap('op', { id: 1 }, expensiveOp);
    await nullCacheWrapper.wrap('op', { id: 1 }, expensiveOp);

    // Both calls should execute
    expect(callCount).toBe(2);
  });

  test('handles null/undefined results', async () => {
    const nullOp = async () => null;
    const undefinedOp = async () => undefined;

    // Should not cache null/undefined
    await wrapper.wrap('op1', { id: 1 }, nullOp);
    await wrapper.wrap('op2', { id: 2 }, undefinedOp);

    const stats = cache.getStats();
    expect(stats.sets).toBe(0);
  });

  test('respects TTL', async () => {
    const expensiveOp = async () => ({ result: 'data' });

    await wrapper.wrap('op', { id: 1 }, expensiveOp, { ttl: 100 });

    // Immediate - cached
    const cached = await cache.get('op:{"id":1}');
    expect(cached).toEqual({ result: 'data' });

    // After TTL - expired
    await new Promise(resolve => setTimeout(resolve, 150));
    const expired = await cache.get('op:{"id":1}');
    expect(expired).toBeNull();
  });

  test('wrapWithKey uses custom cache key', async () => {
    let callCount = 0;
    const expensiveOp = async () => {
      callCount++;
      return { result: 'data' };
    };

    const customKey = 'my-custom-key';

    await wrapper.wrapWithKey(customKey, expensiveOp);
    await wrapper.wrapWithKey(customKey, expensiveOp);

    // Should only call once
    expect(callCount).toBe(1);

    // Verify stored with custom key
    const cached = await cache.get(customKey);
    expect(cached).toEqual({ result: 'data' });
  });
});

describe('CacheWrapper.hashString', () => {
  test('generates consistent hashes', () => {
    const str = 'test string';
    const hash1 = CacheWrapper.hashString(str);
    const hash2 = CacheWrapper.hashString(str);
    expect(hash1).toBe(hash2);
  });

  test('generates different hashes for different strings', () => {
    const hash1 = CacheWrapper.hashString('string1');
    const hash2 = CacheWrapper.hashString('string2');
    expect(hash1).not.toBe(hash2);
  });

  test('handles large strings', () => {
    const largeString = 'x'.repeat(10000);
    const hash = CacheWrapper.hashString(largeString);
    expect(typeof hash).toBe('string');
    expect(hash.length).toBeGreaterThan(0);
  });
});
