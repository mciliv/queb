/**
 * Custom Jest matchers for more expressive tests
 */

const matchers = {
  /**
   * Check if value is a valid SMILES string
   */
  toBeValidSmiles(received) {
    const smilesRegex = /^[A-Za-z0-9@+\-\[\]()=#$\\.\\/:]+$/;
    const pass = typeof received === 'string' && 
                 received.length > 0 && 
                 smilesRegex.test(received);
    
    return {
      pass,
      message: () => pass
        ? `expected "${received}" not to be a valid SMILES string`
        : `expected "${received}" to be a valid SMILES string`
    };
  },

  /**
   * Check if object is valid molecular data
   */
  toBeMolecularData(received) {
    const isValid = received &&
      typeof received === 'object' &&
      typeof received.name === 'string' &&
      typeof received.smiles === 'string' &&
      received.name.length > 0 &&
      this.toBeValidSmiles(received.smiles).pass;
    
    const missing = [];
    if (!received) {
      missing.push('object is null/undefined');
    } else {
      if (!received.name) missing.push('name');
      if (!received.smiles) missing.push('smiles');
      if (received.smiles && !this.toBeValidSmiles(received.smiles).pass) {
        missing.push('valid SMILES format');
      }
    }
    
    return {
      pass: isValid,
      message: () => isValid
        ? `expected object not to be valid molecular data`
        : `expected object to be valid molecular data${missing.length ? ` (missing: ${missing.join(', ')})` : ''}`
    };
  },

  /**
   * Check if response is successful
   */
  toBeSuccessResponse(received) {
    const isSuccess = received &&
      typeof received === 'object' &&
      (received.success === true || 
       received.status === 'success' || 
       received.ok === true ||
       (typeof received.status === 'number' && received.status >= 200 && received.status < 300));
    
    return {
      pass: isSuccess,
      message: () => isSuccess
        ? `expected response not to be successful`
        : `expected response to be successful (got: ${JSON.stringify(received)})`
    };
  },

  /**
   * Check if response is an error
   */
  toBeErrorResponse(received, expectedStatus) {
    const isError = received &&
      typeof received === 'object' &&
      (received.success === false || 
       received.status === 'error' || 
       received.ok === false ||
       (typeof received.status === 'number' && received.status >= 400));
    
    const statusMatches = expectedStatus === undefined || 
                         received.status === expectedStatus;
    
    const pass = isError && statusMatches;
    
    return {
      pass,
      message: () => {
        if (pass) {
          return `expected response not to be an error${expectedStatus ? ` with status ${expectedStatus}` : ''}`;
        }
        if (!isError) {
          return `expected response to be an error (got: ${JSON.stringify(received)})`;
        }
        return `expected error status ${expectedStatus} but got ${received.status}`;
      }
    };
  },

  /**
   * Check if array contains molecular data
   */
  toContainMolecules(received, expectedCount) {
    const isArray = Array.isArray(received);
    const allValid = isArray && received.every(item => 
      this.toBeMolecularData(item).pass
    );
    const countMatches = expectedCount === undefined || 
                        received.length === expectedCount;
    
    const pass = isArray && allValid && countMatches;
    
    return {
      pass,
      message: () => {
        if (!isArray) return `expected an array but got ${typeof received}`;
        if (!allValid) return `expected all items to be valid molecular data`;
        if (!countMatches) return `expected ${expectedCount} molecules but got ${received.length}`;
        return `expected not to contain ${expectedCount || ''} valid molecules`;
      }
    };
  },

  /**
   * Check if value is within range
   */
  toBeWithinRange(received, min, max) {
    const pass = typeof received === 'number' && 
                 received >= min && 
                 received <= max;
    
    return {
      pass,
      message: () => pass
        ? `expected ${received} not to be within range [${min}, ${max}]`
        : `expected ${received} to be within range [${min}, ${max}]`
    };
  },

  /**
   * Check if duration is acceptable
   */
  toCompleteWithin(received, maxMs) {
    if (typeof received !== 'function' && !(received instanceof Promise)) {
      return {
        pass: false,
        message: () => `expected a function or promise but got ${typeof received}`
      };
    }
    
    return {
      pass: false,
      message: () => `Use expect(async () => await someFunction()).toCompleteWithin(ms)`,
      // This is an async matcher
      async: true,
      asyncMatch: async () => {
        const start = Date.now();
        try {
          if (typeof received === 'function') {
            await received();
          } else {
            await received;
          }
          const duration = Date.now() - start;
          const pass = duration <= maxMs;
          
          return {
            pass,
            message: () => pass
              ? `expected to take more than ${maxMs}ms but completed in ${duration}ms`
              : `expected to complete within ${maxMs}ms but took ${duration}ms`
          };
        } catch (error) {
          return {
            pass: false,
            message: () => `expected to complete but threw: ${error.message}`
          };
        }
      }
    };
  },

  /**
   * Check if object matches shape
   */
  toMatchShape(received, shape) {
    const checkShape = (obj, shapeObj, path = '') => {
      const errors = [];
      
      for (const [key, expectedType] of Object.entries(shapeObj)) {
        const currentPath = path ? `${path}.${key}` : key;
        
        if (!(key in obj)) {
          errors.push(`missing property "${currentPath}"`);
          continue;
        }
        
        const value = obj[key];
        const actualType = Array.isArray(value) ? 'array' : typeof value;
        
        if (typeof expectedType === 'string') {
          if (actualType !== expectedType) {
            errors.push(`"${currentPath}" should be ${expectedType} but is ${actualType}`);
          }
        } else if (typeof expectedType === 'object') {
          if (actualType !== 'object') {
            errors.push(`"${currentPath}" should be object but is ${actualType}`);
          } else {
            errors.push(...checkShape(value, expectedType, currentPath));
          }
        }
      }
      
      return errors;
    };
    
    const errors = checkShape(received, shape);
    const pass = errors.length === 0;
    
    return {
      pass,
      message: () => pass
        ? `expected object not to match shape`
        : `expected object to match shape:\n${errors.map(e => `  - ${e}`).join('\n')}`
    };
  },

  /**
   * Check if file exists
   */
  toExistAsFile(received) {
    const fs = require('fs');
    const exists = typeof received === 'string' && fs.existsSync(received);
    const isFile = exists && fs.statSync(received).isFile();
    
    return {
      pass: isFile,
      message: () => {
        if (!exists) return `expected file "${received}" to exist`;
        if (!isFile) return `expected "${received}" to be a file but it's a directory`;
        return `expected file "${received}" not to exist`;
      }
    };
  },

  /**
   * Check if API endpoint is reachable
   */
  async toBeReachable(url, options = {}) {
    try {
      const response = await fetch(url, {
        method: options.method || 'GET',
        timeout: options.timeout || 5000,
        ...options
      });
      
      const pass = response.ok || options.allowedStatuses?.includes(response.status);
      
      return {
        pass,
        message: () => pass
          ? `expected ${url} not to be reachable`
          : `expected ${url} to be reachable but got status ${response.status}`
      };
    } catch (error) {
      return {
        pass: false,
        message: () => `expected ${url} to be reachable but got error: ${error.message}`
      };
    }
  }
};

// Export function to extend Jest with these matchers
const extendExpect = () => {
  expect.extend(matchers);
};

module.exports = {
  matchers,
  extendExpect
};
