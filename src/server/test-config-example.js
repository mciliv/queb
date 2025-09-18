// Example test configuration for Structuralizer
// This allows you to specify exact model/prompt combinations for testing

const testConfigs = {
  // Test with specific model and custom prompt
  customPrompt: {
    model: 'gpt-4o',
    prompt: 'custom' // Uses the original requestParams.messages
  },
  
  // Test with specific model and hardcoded prompt
  hardcodedPrompt: {
    model: 'gpt-4-turbo',
    prompt: 'Analyze this chemical structure and return JSON with chemicals array containing name and SMILES.'
  },
  
  // Test with specific model and structuralize prompt
  structuralizeTest: {
    model: 'gpt-4o',
    prompt: 'Identify all chemical compounds in this text and return JSON with chemicals array.'
  },
  
  // Test with specific model and object detection prompt
  objectDetectionTest: {
    model: 'gpt-4o',
    prompt: 'Detect objects in this image and return JSON with object name and recommended bounding box.'
  }
};

// Usage example:
// const structuralizer = new Structuralizer(apiKey, testConfigs.customPrompt);
// const structuralizer = new Structuralizer(apiKey, testConfigs.hardcodedPrompt);

module.exports = testConfigs;
