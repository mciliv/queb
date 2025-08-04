# Molecular Accuracy Testing

This directory contains comprehensive tests for validating the accuracy of our AtomPredictor prompt engineering improvements against known chemical compositions.

## Overview

We've created a sophisticated testing system that validates:
- **SMILES accuracy**: Correct molecular notation generation
- **Chemical composition**: Realistic molecular breakdowns  
- **Prompt effectiveness**: How well our git history improvements work
- **Error handling**: Graceful fallback behavior

## Test Architecture

### 📁 Test Files

1. **`test/unit/prompt-accuracy.test.js`**
   - Mock-based unit tests for prompt behavior
   - Fast execution (uses mocked OpenAI responses)
   - Tests prompt structure and fallback logic
   - Validates SMILES format and validation

2. **`test/integration/molecular-accuracy.test.js`**
   - Real OpenAI API integration tests
   - Tests actual AtomPredictor responses
   - Requires OPENAI_API_KEY environment variable
   - Comprehensive accuracy scoring system

3. **`test/run-accuracy-tests.js`**
   - Convenient test runner with options
   - Handles API key management
   - Supports filtering by material or test suite

### 🎯 Reference Materials

Our tests validate against scientifically accurate compositions:

```javascript
// Simple compounds - must be 100% accurate
"water": {
  expectedSMILES: ["O"],
  forbiddenSMILES: ["H2O", "HOH"], // Should use SMILES not formulas
  accuracy_threshold: 95
}

// Complex beverages - realistic compositions
"red wine": {
  expectedSMILES: ["CCO", "O", "OC(C(O)C(O)=O)C(O)=O"], // Ethanol, water, tartaric acid
  requiredNames: ["ethanol", "water", "tartaric"],
  minComponents: 3,
  accuracy_threshold: 75
}
```

## Running Tests

### Quick Start

```bash
# Run all accuracy tests (unit + integration if API key available)
npm run test:accuracy

# Run only unit tests (no API key needed)
npm run test:accuracy:unit

# Run only integration tests (needs API key)
OPENAI_API_KEY=sk-your-key npm run test:accuracy:integration
```

### Detailed Options

```bash
# Test specific material
node test/run-accuracy-tests.js --material=water --api-key=sk-xxx

# Verbose output with detailed results
node test/run-accuracy-tests.js --verbose

# Test different suites
node test/run-accuracy-tests.js --suite=unit
node test/run-accuracy-tests.js --suite=integration
node test/run-accuracy-tests.js --suite=all
```

### Standard Jest Commands

```bash
# Run specific test files
npm run test:prompt      # Unit tests only
npm run test:molecular   # Integration tests only

# Watch mode for development
npm run test:watch

# All tests with verbose output
npm run test:verbose
```

## Accuracy Scoring System

Our tests use a comprehensive scoring system (0-100%):

### Score Breakdown
- **Expected SMILES (40 points)**: Contains chemically accurate SMILES
- **Required Names (25 points)**: Includes expected chemical names  
- **No Forbidden SMILES (20 points)**: Avoids chemical formulas like "H2O"
- **Component Count (10 points)**: Appropriate number of chemicals
- **SMILES Validity (5 points)**: Proper SMILES syntax

### Accuracy Thresholds
- **Simple compounds**: 95% (water, ethanol, salt)
- **Beverages**: 75-80% (wine, coffee, beer)
- **Materials**: 70-75% (plastics, minerals)
- **Biological**: 65% (fruits, organic matter)

## Understanding Test Output

### Example Output
```
=== WATER ===
Score: 95% (threshold: 95%)
Breakdown: {
  expectedSMILES: "✓ Found",
  requiredNames: "1/1 found", 
  forbiddenSMILES: "✓ None found",
  componentCount: "1 components",
  smilesValidity: "✓ All valid"
}
Chemicals found: Water: O
```

### Performance Metrics
```
=== PROMPT ENGINEERING EFFECTIVENESS ===
Average accuracy: 78.2%
Pass rate: 80.0%
Average chemicals per analysis: 3.4
Materials tested: 10
```

## What We're Testing

### 🔬 Chemical Accuracy
- **Correct SMILES notation**: "O" not "H2O"
- **Ionic compounds**: "[Na+].[Cl-]" not "NaCl"
- **Complex molecules**: Proper representation with length limits
- **Realistic compositions**: Wine has ethanol + water + acids

### 🎯 Prompt Engineering Effectiveness
- **Git history improvements**: Best techniques from working versions
- **Context awareness**: Different examples for beverages vs materials
- **Fallback quality**: Smart error handling for edge cases
- **Validation**: SMILES format checking and quality control

### 🚫 Common Error Prevention
- **Formula errors**: Chemical formulas instead of SMILES
- **Overly complex SMILES**: Length and complexity constraints
- **Missing components**: Unrealistic single-chemical responses
- **Invalid syntax**: Malformed SMILES strings

## Integration with Development

### Continuous Testing
```bash
# Before committing prompt changes
npm run test:accuracy:unit

# Before deploying (with API key)
npm run test:accuracy
```

### Monitoring Improvements
The tests track accuracy metrics over time, helping validate that prompt engineering changes actually improve results.

### Adding New Materials
```javascript
// In test/integration/molecular-accuracy.test.js
"new_material": {
  expectedSMILES: ["expected_smiles_here"],
  requiredNames: ["chemical_names"],
  accuracy_threshold: 70
}
```

## Troubleshooting

### No API Key Issues
```
⚠️  Skipping molecular accuracy tests - no OPENAI_API_KEY set
```
**Solution**: Set `OPENAI_API_KEY` environment variable or use `--api-key` option

### Low Accuracy Scores
```
✗ water: 45% (need 95%)
Issues: Missing expected SMILES: O
```
**Solution**: Check prompt instructions, examples, or constraints in `backend/prompts/`

### Test Timeouts
```
Timeout - Async callback was not invoked within the 30000ms timeout
```
**Solution**: API calls can be slow. Tests have 30s timeouts for single materials, 2min for full suites.

## Development Guidelines

### When to Run Tests
- ✅ After modifying prompt instructions
- ✅ Before committing AtomPredictor changes
- ✅ When adding new chemical examples
- ✅ Before production deployments

### Interpreting Results
- **>90% accuracy**: Excellent, ready for production
- **70-90% accuracy**: Good, monitor for improvements  
- **50-70% accuracy**: Acceptable for complex materials
- **<50% accuracy**: Needs prompt engineering work

### Adding Test Cases
1. Research actual chemical composition
2. Define expected SMILES and names
3. Set appropriate accuracy threshold
4. Test with both unit and integration suites

## Scientific Accuracy

Our test materials are based on:
- **PubChem database**: Verified SMILES notation
- **Scientific literature**: Realistic compositions
- **Chemistry databases**: Standard molecular representations
- **Git history analysis**: Previously working examples

This ensures tests validate real chemical accuracy, not just format compliance.

---

## Quick Reference

| Command | Purpose |
|---------|---------|
| `npm run test:accuracy` | Full accuracy test suite |
| `npm run test:accuracy:unit` | Fast unit tests only |
| `npm run test:accuracy:integration` | Real API tests only |
| `node test/run-accuracy-tests.js --help` | Show all options |
| `npm run test:watch` | Development watch mode |

**Need help?** Check the test output for specific accuracy breakdowns and issues. 