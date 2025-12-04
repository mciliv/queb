# Backend Text Validation Tests

This test file (`structuralize-text-validation.test.js`) contains comprehensive tests for the `/api/structuralize` endpoint's text parameter validation.

## Running the Tests

### Option 1: Run with Jest directly
```bash
npx jest tests/integration/structuralize-text-validation.test.js --config tests/jest.integration.js
```

### Option 2: Run all integration tests
```bash
npm run test:integration
```

### Option 3: Run with specific test pattern
```bash
npm run test:integration -- --testNamePattern="Text Parameter Validation"
```

## Test Coverage

The tests verify:

### ✅ Valid Inputs
- Valid text strings
- Text with whitespace (should be trimmed)
- Default lookupMode handling

### ❌ Invalid Inputs (should return 400)
- Missing text parameter
- Empty string
- Whitespace-only string
- Null value
- Undefined value
- Non-string types (number, boolean, object, array)
- Empty request body

### 🔍 Edge Cases
- Very long text strings
- Special characters
- Unicode characters

### 📋 lookupMode Validation
- Valid lookupMode
- Invalid lookupMode type

## Expected Behavior

All invalid inputs should:
1. Return HTTP 400 status
2. Include error message: "Missing or invalid 'text' parameter"
3. NOT call the OpenAI API (mocked service should not be invoked)

Valid inputs should:
1. Return HTTP 200 status
2. Process the request through the full pipeline
3. Call the OpenAI API with trimmed text
