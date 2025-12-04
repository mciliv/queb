// Google Cloud Secret Manager utility
const { SecretManagerServiceClient } = require('@google-cloud/secret-manager');

class SecretManager {
  constructor() {
    this.client = new SecretManagerServiceClient();
    this.projectId = process.env.GOOGLE_CLOUD_PROJECT || 'mol-analysis-app';
    this.cache = new Map();
    this.cacheExpiry = 5 * 60 * 1000; // 5 minutes
  }

  /**
   * Get secret value from Google Cloud Secret Manager
   * @param {string} secretName - Name of the secret
   * @param {string} version - Version of the secret (default: 'latest')
   * @returns {Promise<string>} Secret value
   */
  async getSecret(secretName, version = 'latest') {
    const cacheKey = `${secretName}:${version}`;
    
    // Check cache first
    if (this.cache.has(cacheKey)) {
      const cached = this.cache.get(cacheKey);
      if (Date.now() - cached.timestamp < this.cacheExpiry) {
        return cached.value;
      }
      this.cache.delete(cacheKey);
    }

    try {
      const name = `projects/${this.projectId}/secrets/${secretName}/versions/${version}`;
      const [response] = await this.client.accessSecretVersion({ name });
      const secretValue = response.payload.data.toString();
      
      // Cache the result
      this.cache.set(cacheKey, {
        value: secretValue,
        timestamp: Date.now()
      });
      
      return secretValue;
    } catch (error) {
      console.error(`Error accessing secret ${secretName}:`, error.message);
      throw new Error(`Failed to access secret: ${secretName}`);
    }
  }

  /**
   * Get OpenAI API key based on environment
   * Production uses 'OPENAI' secret, Development uses 'openai-api-key' secret
   * @returns {Promise<string>} OpenAI API key
   */
  async getOpenAIKey() {
    const isProduction = process.env.NODE_ENV === 'production';
    const secretName = isProduction ? 'OPENAI' : 'openai-api-key';
    return await this.getSecret(secretName);
  }

  /**
   * Clear the cache (useful for testing or forced refresh)
   */
  clearCache() {
    this.cache.clear();
  }
}

// Export singleton instance
const secretManager = new SecretManager();
module.exports = secretManager;

