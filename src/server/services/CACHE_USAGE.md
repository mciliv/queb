# Cache System Usage Guide

## Overview

Extracted from `Structuralizer.js` (lines 102-119, 291-317) per user request:
> "extract the cache components into it's own file, using the simplest way to integrate with existing code when needed"

**Components:**
- `SimpleCache.js` - In-memory cache with TTL support
- `CacheWrapper.js` - Function wrapper for caching expensive operations

**Why not use a library?**
- No external dependencies needed (no Redis, node-cache, etc.)
- ~200 lines of simple, maintainable code
- Tailored to your exact needs (payload hashing, flexible TTL)
- Can be replaced with Redis later without changing consumer code

## Quick Start

### 1. Basic Setup

```javascript
const SimpleCache = require('./services/SimpleCache');
const CacheWrapper = require('./services/CacheWrapper');

// Create cache instance
const cache = new SimpleCache({
  maxSize: 1000,        // Max number of entries
  defaultTTL: 300000,   // 5 minutes
  cleanupInterval: 60000 // Cleanup every minute
});

// Create wrapper
const cacheWrapper = new CacheWrapper({
  cache: cache,
  enabled: true,
  logger: console,
  defaultTTL: 300000
});
```

### 2. Wrap Expensive Operations

```javascript
// Before: Direct call
const result = await expensiveOperation(params);

// After: With caching
const result = await cacheWrapper.wrap(
  'operation-name',      // Cache key prefix
  { userId: 123 },       // Data for cache key
  () => expensiveOperation(params),  // Function to cache
  { ttl: 300000 }        // Optional TTL override
);
```

## Integration with Structuralizer

### Current Code (Structuralizer.js lines 102-119)

```javascript
async chemicals(payload) {
  const startTime = Date.now();

  if (this.cache && this.config.cacheEnabled) {
    const cacheKey = this._getCacheKey(payload);
    const cached = await this.cache.get(cacheKey);
    if (cached) {
      this.logger.info('Cache hit for prediction', {
        object: payload.object,
        duration: Date.now() - startTime
      });
      return cached;
    }
  }

  const result = await this._chemicals(payload.object);

  if (this.cache && this.config.cacheEnabled && result) {
    const cacheKey = this._getCacheKey(payload);
    await this.cache.set(cacheKey, result, 300000);
  }

  return result;
}
```

### Updated Code (Using CacheWrapper)

```javascript
const CacheWrapper = require('./CacheWrapper');

class Structuralizer {
  constructor(dependencies = {}) {
    // ... existing code ...
    this.cacheWrapper = new CacheWrapper({
      cache: dependencies.cache,
      enabled: dependencies.config?.cacheEnabled,
      logger: dependencies.logger
    });
  }

  async chemicals(payload) {
    return await this.cacheWrapper.wrap(
      'structuralizer',
      {
        object: payload.object,
        lookupMode: payload.lookupMode,
        x: payload.x,
        y: payload.y,
        imageHash: payload.imageBase64 
          ? CacheWrapper.hashString(payload.imageBase64.substring(0, 100))
          : null
      },
      () => this._chemicals(payload),
      { ttl: 300000 }
    );
  }
  
  // Remove _getCacheKey and _hashString methods - now in CacheWrapper
}
```

### Update server.js

```javascript
const SimpleCache = require('./src/server/services/SimpleCache');

// Create cache instance
const cache = new SimpleCache({
  maxSize: 1000,
  defaultTTL: 300000
});

const structuralizer = new Structuralizer({
  aiService,
  molecularProcessor,
  nameResolver: { resolveName },
  promptEngine,
  errorHandler,
  logger,
  cache: cache,  // <-- Add cache
  config: {
    aiTimeout: 10000,
    maxRetries: 2,
    cacheEnabled: true  // <-- Enable caching
  }
});
```

## Advanced Usage

### Caching with Large Payloads (Images)

```javascript
const CacheWrapper = require('./CacheWrapper');

async function analyzeImage(imageBase64, metadata) {
  return await cacheWrapper.wrap(
    'image-analysis',
    {
      ...metadata,
      // Hash first 100 chars instead of storing full base64 in key
      imageHash: CacheWrapper.hashString(imageBase64.substring(0, 100))
    },
    () => expensiveImageAnalysis(imageBase64, metadata),
    { ttl: 600000 } // 10 minutes
  );
}
```

### Multiple Cache Instances

```javascript
// Short-term cache for API responses
const apiCache = new SimpleCache({
  maxSize: 500,
  defaultTTL: 60000  // 1 minute
});

// Long-term cache for expensive computations
const computeCache = new SimpleCache({
  maxSize: 100,
  defaultTTL: 3600000  // 1 hour
});

const apiWrapper = new CacheWrapper({ cache: apiCache });
const computeWrapper = new CacheWrapper({ cache: computeCache });
```

### Monitoring Cache Performance

```javascript
// Get cache statistics
const stats = cache.getStats();
console.log(stats);
// {
//   hits: 150,
//   misses: 50,
//   sets: 50,
//   evictions: 0,
//   size: 50,
//   hitRate: '75.00%'
// }
```

### Graceful Shutdown

```javascript
process.on('SIGTERM', async () => {
  cache.destroy();  // Stop cleanup timer
  await cache.clear();  // Clear all entries
  process.exit(0);
});
```

## Decorator Pattern (Not Recommended)

User asked about decorators. Here's why we didn't use them:

**Option 1: Babel Decorators (requires build step)**
```javascript
// Requires @babel/plugin-proposal-decorators
class MyService {
  @cached('my-operation', { ttl: 60000 })
  async expensiveMethod(params) {
    // ...
  }
}
```

**Problems:**
- JavaScript decorators are still Stage 3 proposal (not stable)
- Requires Babel/TypeScript transpilation
- Adds build complexity
- Less clear for developers unfamiliar with decorators
- Harder to debug

**Our Solution: Simple Function Wrapper**
```javascript
async expensiveMethod(params) {
  return await this.cacheWrapper.wrap(
    'my-operation',
    params,
    () => this._expensiveMethod(params),
    { ttl: 60000 }
  );
}
```

**Benefits:**
- Works in any JavaScript environment (no build tools)
- Clear and explicit
- Easy to debug
- Flexible (can customize per call)

## Standalone Package

This cache system can be extracted into a standalone npm package:

**Package Structure:**
```
simple-cache-wrapper/
├── package.json
├── README.md
├── src/
│   ├── SimpleCache.js
│   └── CacheWrapper.js
└── tests/
    ├── SimpleCache.test.js
    └── CacheWrapper.test.js
```

**To make it standalone:**

1. **Add TypeScript definitions** (simple-cache-wrapper.d.ts)
2. **Add comprehensive tests**
3. **Document cache interface contract:**
   ```typescript
   interface Cache {
     get(key: string): Promise<any | null>;
     set(key: string, value: any, ttl: number): Promise<void>;
     delete(key: string): Promise<boolean>;
     clear(): Promise<void>;
   }
   ```
4. **Publish to npm** as `@queb/simple-cache-wrapper`

**Benefits of keeping it in queb:**
- No extra dependency to maintain
- Can iterate quickly
- ~200 lines is manageable

**When to extract:**
- When other projects need it
- When it needs significant features (Redis adapter, etc.)
- When you want community contributions

## Testing

```javascript
// tests/unit/cache-wrapper.test.js
const SimpleCache = require('../../src/server/services/SimpleCache');
const CacheWrapper = require('../../src/server/services/CacheWrapper');

describe('CacheWrapper', () => {
  test('caches expensive operations', async () => {
    const cache = new SimpleCache();
    const wrapper = new CacheWrapper({ cache });
    
    let callCount = 0;
    const expensiveOp = async () => {
      callCount++;
      return { result: 'data' };
    };
    
    // First call - cache miss
    const result1 = await wrapper.wrap('test', { id: 1 }, expensiveOp);
    expect(callCount).toBe(1);
    
    // Second call - cache hit
    const result2 = await wrapper.wrap('test', { id: 1 }, expensiveOp);
    expect(callCount).toBe(1); // Not called again
    expect(result2).toEqual(result1);
  });
  
  test('respects TTL', async () => {
    const cache = new SimpleCache();
    const wrapper = new CacheWrapper({ cache });
    
    await wrapper.wrap('test', { id: 1 }, async () => 'data', { ttl: 100 });
    
    // Immediate read - cached
    const cached = await cache.get('test:{"id":1}');
    expect(cached).toBe('data');
    
    // After TTL - expired
    await new Promise(resolve => setTimeout(resolve, 150));
    const expired = await cache.get('test:{"id":1}');
    expect(expired).toBeNull();
  });
});
```

## Migration Path

**Current State:** Cache disabled, no infrastructure
**Step 1:** Enable SimpleCache (this PR)
**Step 2:** Monitor hit rates, adjust TTLs
**Step 3:** (Optional) Migrate to Redis if needed:

```javascript
// Drop-in replacement with Redis
const redis = require('redis');
const client = redis.createClient();

const redisCache = {
  async get(key) {
    return JSON.parse(await client.get(key));
  },
  async set(key, value, ttl) {
    await client.setEx(key, Math.floor(ttl / 1000), JSON.stringify(value));
  },
  async delete(key) {
    await client.del(key);
  },
  async clear() {
    await client.flushDb();
  }
};

// Same CacheWrapper code works!
const wrapper = new CacheWrapper({ cache: redisCache });
```

No changes to Structuralizer or other consumers needed.
