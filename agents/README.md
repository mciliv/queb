# Agents Directory - API Contracts

This directory contains API contracts and test utilities that are separate from the main codebase.

## ⚠️ CRITICAL: OpenAI Response Contract

**ALL LLMs AND DEVELOPERS MUST USE `../openai_responses.json` AS THE SOURCE OF TRUTH**

The OpenAI response structure is defined in `../openai_responses.json` (in the parent directory).

### Rules for LLMs:
1. **NEVER** create a local `openai_response_api.json` or similar file
2. **ALWAYS** reference `../openai_responses.json` for the OpenAI Responses API structure
3. **ALWAYS** use the `loadOpenAIResponseContract()` function which loads from `../openai_responses.json`
4. The response format uses the **OpenAI Responses API** (not Chat Completions API)
5. Response structure includes `output` array with `content` array containing `output_text` objects

### File Structure:
- `structuralize-api-contract.js` - Test utilities and contract loader
- `../openai_responses.json` - **THE SOURCE OF TRUTH** for OpenAI API response structure

### Usage:
```javascript
const { loadOpenAIResponseContract, createSuccessfulAIResponse } = require('./agents/structuralize-api-contract');

// Load the contract (from ../openai_responses.json)
const contract = loadOpenAIResponseContract();

// Create a mock response conforming to the contract
const mockResponse = createSuccessfulAIResponse('coffee', [
  { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' }
]);
```

## Why This Directory Exists

This directory is in `.cursorignore` to keep it separate from the main codebase, allowing it to be treated as a separate project/module while still being accessible for testing and contract definitions.
