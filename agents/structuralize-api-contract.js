/**
 * API Contract for /api/structuralize Endpoint
 * 
 * This file defines the contract and test setup utilities for testing
 * the /api/structuralize endpoint. It provides a standardized way to
 * create test containers and mocks for consistent testing across the project.
 * 
 * This file is in agents/ to keep it separate from the main codebase
 * and allow it to be treated as a separate project/module.
 * 
 * ⚠️ CRITICAL: The OpenAI response structure MUST be loaded from ../openai_responses.json
 * This is the source of truth for the OpenAI Responses API format.
 * DO NOT create or use a local openai_response_api.json file.
 * ALWAYS reference ../openai_responses.json for the response contract.
 */

const fs = require('fs');
const path = require('path');

/**
 * Loads the OpenAI response API contract from JSON
 * 
 * ⚠️ IMPORTANT: This function MUST load from ../openai_responses.json
 * This is the authoritative source for OpenAI Responses API structure.
 * LLMs and developers MUST use this file, not any local copy.
 * 
 * @returns {Object} The contract schema from ../openai_responses.json
 */
function loadOpenAIResponseContract() {
  // CRITICAL: Use ../openai_responses.json - this is the source of truth
  const contractPath = path.join(__dirname, '..', '..', 'openai_responses.json');
  
  if (!fs.existsSync(contractPath)) {
    throw new Error(
      `OpenAI response contract not found at ${contractPath}. ` +
      `This file MUST exist and MUST be used as the source of truth for OpenAI Responses API structure.`
    );
  }
  
  try {
    const contract = JSON.parse(fs.readFileSync(contractPath, 'utf8'));
    return contract;
  } catch (error) {
    throw new Error(`Failed to load OpenAI response contract from ${contractPath}: ${error.message}`);
  }
}

// Load the contract - this MUST come from ../openai_responses.json
const OPENAI_RESPONSE_CONTRACT = loadOpenAIResponseContract();

/**
 * Creates a mock OpenAI client structure for testing
 * @returns {Object} Mock OpenAI client with chat.completions.create method
 */
function createOpenAIMock() {
  return {
    openaiClient: {
      chat: {
        completions: {
          create: jest.fn()
        }
      }
    }
  };
}

/**
 * Creates a successful AI response mock based on the contract from ../openai_responses.json
 * 
 * ⚠️ IMPORTANT: This function MUST conform to the structure in ../openai_responses.json
 * The response format uses the OpenAI Responses API (not Chat Completions API)
 * 
 * @param {string} object - The object name
 * @param {Array} chemicals - Array of chemical objects with name and smiles
 * @param {Object} options - Additional options for the response
 * @returns {Object} Mock OpenAI API response conforming to ../openai_responses.json structure
 */
function createSuccessfulAIResponse(object, chemicals = [], options = {}) {
  const contract = OPENAI_RESPONSE_CONTRACT;
  const parsedContent = {
    object,
    chemicals
  };

  // Use the Responses API format from ../openai_responses.json
  // Structure: output array with content array containing output_text objects
  return {
    id: options.id || `resp_${Date.now()}`,
    object: contract.object || 'response',
    created_at: options.created_at || Math.floor(Date.now() / 1000),
    status: options.status || 'completed',
    error: null,
    incomplete_details: null,
    instructions: null,
    max_output_tokens: null,
    model: options.model || contract.model || 'gpt-4.1-2025-04-14',
    output: [
      {
        type: 'message',
        id: options.messageId || `msg_${Date.now()}`,
        status: 'completed',
        role: options.role || 'assistant',
        content: [
          {
            type: 'output_text',
            text: JSON.stringify(parsedContent),
            annotations: []
          }
        ]
      }
    ],
    parallel_tool_calls: options.parallel_tool_calls !== undefined ? options.parallel_tool_calls : true,
    previous_response_id: null,
    reasoning: {
      effort: null,
      summary: null
    },
    store: options.store !== undefined ? options.store : true,
    temperature: options.temperature || contract.temperature || 1.0,
    text: {
      format: {
        type: 'text'
      }
    },
    tool_choice: options.tool_choice || 'auto',
    tools: options.tools || [],
    top_p: options.top_p || contract.top_p || 1.0,
    truncation: options.truncation || 'disabled',
    usage: options.usage || {
      input_tokens: 100,
      input_tokens_details: {
        cached_tokens: 0
      },
      output_tokens: 50,
      output_tokens_details: {
        reasoning_tokens: 0
      },
      total_tokens: 150
    },
    user: null,
    metadata: options.metadata || {}
  };
}

/**
 * Test setup contract for /api/structuralize endpoint tests
 * 
 * This provides a standardized way to set up tests with:
 * - Mocked OpenAI client
 * - Test container
 * - Express app instance
 * 
 * @param {Function} createTestContainer - Function to create test container
 * @param {Function} createApp - Function to create Express app
 * @returns {Object} Test setup utilities
 */
function createTestSetup(createTestContainer, createApp) {
  return {
    /**
     * Sets up test environment before each test
     * @returns {Promise<Object>} Object with app, container, and mocks
     */
    async setup() {
      const mocks = createOpenAIMock();
      const container = createTestContainer(mocks);
      const app = await createApp(container);
      
      return { app, container, mocks };
    },

    /**
     * Cleans up test environment after each test
     * @param {Object} testContext - Object with container from setup()
     */
    teardown(testContext) {
      if (testContext && testContext.container) {
        testContext.container.clear();
      }
      jest.clearAllMocks();
    },

    /**
     * Creates a mock for a successful structuralize response
     * @param {string} objectName - Name of the object being analyzed
     * @param {Array} chemicals - Array of chemical objects
     */
    mockSuccessfulResponse(objectName, chemicals) {
      return createSuccessfulAIResponse(objectName, chemicals);
    }
  };
}

/**
 * Expected error messages for validation failures
 */
const ERROR_MESSAGES = {
  MISSING_TEXT: 'Missing or invalid "text" parameter',
  INVALID_LOOKUP_MODE: 'Invalid "lookupMode" parameter'
};

/**
 * HTTP status codes used by the API
 */
const STATUS_CODES = {
  SUCCESS: 200,
  BAD_REQUEST: 400,
  SERVER_ERROR: 500
};

module.exports = {
  createOpenAIMock,
  createSuccessfulAIResponse,
  createTestSetup,
  loadOpenAIResponseContract,
  OPENAI_RESPONSE_CONTRACT,
  ERROR_MESSAGES,
  STATUS_CODES
};
