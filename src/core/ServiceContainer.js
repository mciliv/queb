class ServiceContainer {
  constructor() {
    this.services = new Map();
    this.factories = new Map();
    this.resolving = new Set(); // Track resolving services to detect circular deps
  }

  /**
   * Register a service factory
   * @param {string} name - Service identifier
   * @param {Function|*} factory - Factory function or instance
   * @param {Object} options - Registration options
   * @param {boolean} options.singleton - Create only one instance (default: true)
   * @param {Array<string>} options.tags - Service tags for bulk operations
   */
  register(name, factory, options = {}) {
    const config = {
      factory,
      singleton: options.singleton !== false,
      tags: options.tags || [],
      instance: null
    };

    // If factory is not a function, treat it as a value
    if (typeof factory !== 'function') {
      config.instance = factory;
      config.singleton = true;
    }

    this.factories.set(name, config);
    return this;
  }

  /**
   * Register multiple services at once
   * @param {Object} providers - Map of service names to factories
   */
  registerAll(providers) {
    Object.entries(providers).forEach(([name, factory]) => {
      this.register(name, factory);
    });
    return this;
  }

  /**
   * Get a service instance
   * @param {string} name - Service identifier
   * @returns {*} Service instance
   * @throws {Error} If service not found or circular dependency detected
   */
  async get(name) {
    if (!this.factories.has(name)) {
      throw new Error(`Service '${name}' not registered`);
    }

    const config = this.factories.get(name);

    // Return existing singleton instance
    if (config.singleton && config.instance !== null) {
      return config.instance;
    }

    // Check for circular dependencies
    if (this.resolving.has(name)) {
      const cycle = Array.from(this.resolving).join(' -> ') + ' -> ' + name;
      throw new Error(`Circular dependency detected: ${cycle}`);
    }

    // Mark as resolving
    this.resolving.add(name);

    try {
      // Create instance
      let instance;
      if (typeof config.factory === 'function') {
        instance = await config.factory(this);
      } else {
        instance = config.factory;
      }

      // Store singleton
      if (config.singleton) {
        config.instance = instance;
      }

      return instance;
    } finally {
      // Clear resolving flag
      this.resolving.delete(name);
    }
  }

  /**
   * Check if a service is registered
   * @param {string} name - Service identifier
   * @returns {boolean}
   */
  has(name) {
    return this.factories.has(name);
  }

  /**
   * Get all services with a specific tag
   * @param {string} tag - Tag to filter by
   * @returns {Promise<Array>} Array of service instances
   */
  async getTagged(tag) {
    const services = [];
    for (const [name, config] of this.factories) {
      if (config.tags.includes(tag)) {
        services.push(await this.get(name));
      }
    }
    return services;
  }

  /**
   * Create a scoped container (useful for request-scoped services)
   * @returns {ServiceContainer} New container with same factories
   */
  createScope() {
    const scoped = new ServiceContainer();
    // Copy factory definitions but not instances
    for (const [name, config] of this.factories) {
      scoped.factories.set(name, {
        ...config,
        instance: null
      });
    }
    return scoped;
  }

  /**
   * Clear all singleton instances (useful for testing)
   */
  clear() {
    for (const config of this.factories.values()) {
      config.instance = null;
    }
    this.resolving.clear();
  }

  /**
   * Get dependency graph for visualization/debugging
   * @returns {Object} Dependency graph
   */
  async getDependencyGraph() {
    const graph = {};
    
    for (const [name, config] of this.factories) {
      if (typeof config.factory === 'function') {
        // Analyze factory function to find dependencies
        const funcStr = config.factory.toString();
        const containerCalls = funcStr.matchAll(/container\.get\(['"`](\w+)['"`]\)/g);
        const deps = Array.from(containerCalls, m => m[1]);
        graph[name] = deps;
      } else {
        graph[name] = [];
      }
    }
    
    return graph;
  }
}

module.exports = ServiceContainer;
