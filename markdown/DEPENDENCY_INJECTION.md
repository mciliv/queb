# Dependency Injection Implementation Guide

This guide explains how to use the new dependency injection (DI) system in the Queb application, which makes the codebase more testable, maintainable, and flexible.

## Overview

The DI system consists of three main components:

1. **ServiceContainer** - Manages service registration and resolution
2. **ServiceProvider** - Configures all application services
3. **Refactored Services** - Services that accept dependencies through constructors

## Benefits

### 1. **Easier Testing**
```javascript
// Before DI - Complex mocking required
jest.mock('../services/molecular-processor');
jest.mock('../ai/openai/client');
// ... many more mocks

// After DI - Simple mock injection
const mocks = {
  molecularProcessor: { generateSDF: jest.fn() },
  aiClient: { chat: { completions: { create: jest.fn() } } }
};
const service = new Structuralizer(mocks);
```

### 2. **Clear Dependencies**
All service dependencies are explicit in the constructor, making it easy to understand what each service needs.

### 3. **Configuration Flexibility**
Different environments can have different service configurations without changing code.

### 4. **No Hidden State**
Services don't use global imports, eliminating hidden dependencies and singleton issues.

## Usage

### Basic Service Registration

```javascript
const { createContainer } = require('./core/ServiceProvider');

// Create application container
const container = createContainer();

// Get a service
const structuralizer = await container.get('structuralizer');
```

### Creating New Services

When creating a new service, follow this pattern:

```javascript
class MyNewService {
  constructor(dependencies = {}) {
    // Required dependencies
    this.database = dependencies.database;
    this.logger = dependencies.logger;
    
    // Optional dependencies with defaults
    this.cache = dependencies.cache || null;
    
    // Validate required dependencies
    if (!this.database || !this.logger) {
      throw new Error('Missing required dependencies: database, logger');
    }
  }
  
  async doSomething() {
    this.logger.info('Doing something');
    
    // Use cache if available
    if (this.cache) {
      const cached = await this.cache.get('key');
      if (cached) return cached;
    }
    
    // Do actual work
    const result = await this.database.query('...');
    
    if (this.cache) {
      await this.cache.set('key', result);
    }
    
    return result;
  }
}
```

### Registering the Service

Add to `ServiceProvider.js`:

```javascript
container.register('myNewService', async (c) => {
  return new MyNewService({
    database: await c.get('database'),
    logger: await c.get('logger'),
    cache: await c.get('cache')
  });
}, { 
  tags: ['business'],
  singleton: true  // Create only one instance
});
```

## Testing with DI

### Unit Tests

```javascript
describe('MyNewService', () => {
  it('should use cache when available', async () => {
    const mocks = {
      database: { query: jest.fn() },
      logger: { info: jest.fn() },
      cache: { 
        get: jest.fn().mockResolvedValue('cached-result'),
        set: jest.fn()
      }
    };
    
    const service = new MyNewService(mocks);
    const result = await service.doSomething();
    
    expect(result).toBe('cached-result');
    expect(mocks.database.query).not.toHaveBeenCalled();
  });
});
```

### Integration Tests

```javascript
describe('Integration', () => {
  it('should work with real database', async () => {
    // Use test container with only external services mocked
    const container = createTestContainer({
      openaiClient: mockAIClient  // Mock only external service
    });
    
    // Real database and other services will be used
    const app = await createApp(container);
    
    // Test with real integration
    const response = await request(app)
      .post('/api/endpoint')
      .send({ data: 'test' });
      
    expect(response.status).toBe(200);
  });
});
```

## Migration Path

To migrate existing code to use DI:

### 1. Start with Leaf Services
Begin with services that have few dependencies:
- ErrorHandler
- Logger  
- Configuration

### 2. Move Up the Dependency Tree
Next, migrate services that depend on the leaf services:
- NameResolver
- MolecularProcessor
- PromptEngine

### 3. Migrate Business Services
Then migrate higher-level services:
- Structuralizer
- MolecularPredictionService

### 4. Update Entry Points
Finally, update server.js and other entry points to use the container.

### Example Migration

**Before:**
```javascript
// Old service with hidden dependencies
const config = require('../config');
const logger = require('./logger');

class OldService {
  async process(data) {
    logger.info('Processing');
    return config.get('enabled') ? data : null;
  }
}

module.exports = new OldService();  // Singleton
```

**After:**
```javascript
// New service with explicit dependencies
class NewService {
  constructor(dependencies = {}) {
    this.config = dependencies.config;
    this.logger = dependencies.logger;
    
    if (!this.config || !this.logger) {
      throw new Error('Missing required dependencies');
    }
  }
  
  async process(data) {
    this.logger.info('Processing');
    return this.config.get('enabled') ? data : null;
  }
}

module.exports = NewService;  // Export class, not instance
```

## Advanced Features

### Scoped Containers

For request-scoped services:

```javascript
app.use((req, res, next) => {
  // Create scoped container for this request
  req.container = container.createScope();
  
  // Register request-specific services
  req.container.register('requestId', () => req.id);
  req.container.register('user', () => req.user);
  
  next();
});
```

### Service Tags

Group related services:

```javascript
// Get all services tagged as 'external'
const externalServices = await container.getTagged('external');

// Useful for:
// - Bulk initialization
// - Health checks
// - Graceful shutdown
```

### Circular Dependency Detection

The container automatically detects circular dependencies:

```javascript
// This will throw an error
container.register('serviceA', async (c) => {
  const serviceB = await c.get('serviceB');
  return new ServiceA(serviceB);
});

container.register('serviceB', async (c) => {
  const serviceA = await c.get('serviceA');  // Circular!
  return new ServiceB(serviceA);
});
```

## Best Practices

1. **Keep Constructors Simple**: Only assign dependencies, don't do work
2. **Validate Dependencies**: Check required dependencies in constructor
3. **Use Interfaces**: Define clear interfaces for dependencies (consider TypeScript)
4. **Avoid Service Locator**: Don't pass the container itself to services
5. **Prefer Constructor Injection**: Avoid setter injection or property injection
6. **Mock at Boundaries**: Mock external services, use real internal services in tests

## Common Patterns

### Optional Dependencies
```javascript
class Service {
  constructor(deps = {}) {
    this.required = deps.required;
    this.optional = deps.optional || null;
    
    if (!this.required) {
      throw new Error('Required dependency missing');
    }
  }
  
  async doWork() {
    if (this.optional) {
      // Use optional service if available
      await this.optional.enhance();
    }
    return this.required.process();
  }
}
```

### Factory Pattern
```javascript
container.register('serviceFactory', async (c) => {
  const config = await c.get('config');
  
  return {
    create: (type) => {
      switch (type) {
        case 'ai':
          return new AIService({ config });
        case 'mock':
          return new MockService({ config });
        default:
          throw new Error(`Unknown type: ${type}`);
      }
    }
  };
});
```

### Decorator Pattern
```javascript
container.register('cachedStructuralizer', async (c) => {
  const base = await c.get('structuralizer');
  const cache = await c.get('cache');
  
  return new CachedStructuralizer(base, cache);
});
```

## Troubleshooting

### Service Not Found
```
Error: Service 'xyz' not registered
```
- Check spelling in ServiceProvider.js
- Ensure service is registered before use

### Circular Dependency
```
Error: Circular dependency detected: serviceA -> serviceB -> serviceA
```
- Refactor to break the cycle
- Consider using events or callbacks instead

### Missing Dependencies
```
Error: Missing required dependencies: logger, config
```
- Ensure all required services are registered
- Check that dependencies are passed to constructor

## Conclusion

The dependency injection system makes the Queb application more maintainable and testable. By making dependencies explicit and injectable, we can:

- Write faster, more focused tests
- Easily swap implementations
- Understand service relationships
- Avoid hidden coupling
- Scale the application confidently

Start with new services and gradually migrate existing ones. The investment in DI pays off quickly through easier testing and better code clarity.
