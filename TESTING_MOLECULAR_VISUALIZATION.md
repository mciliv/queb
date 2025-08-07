# Molecular Visualization Testing Guide

## Overview
Test the complete molecular visualization pipeline without manual typing using multiple automated approaches.

## Testing Options

### 🎯 1. Full Pipeline Test (Recommended)
**Uses the same chemical data as existing LLM endpoint tests**

```bash
npm run test:pipeline
```

**What it tests:**
- ✅ AI Analysis: Object → Molecules + SMILES
- ✅ SDF Generation: SMILES → 3D structure files  
- ✅ File Serving: SDF files accessible via HTTP
- ✅ Visualization Data: Complete data structure for 3DMol.js
- ✅ Quality Validation: SMILES quality, expected chemicals
- ✅ Error Handling: Invalid inputs, SDF generation failures
- ✅ Performance: Pipeline completion within reasonable time

**Test Categories (Minimal Subset Validation):**
- **Basics**: water, ethanol, sodium chloride (high confidence subset)
- **Beverages**: red wine (ethanol+water), black coffee (water+caffeine) (essential subset)  
- **Biological**: fresh apple (water only) (most confident subset)

**Requirements:**
- Set `OPENAI_API_KEY` environment variable
- Dev server NOT required (creates own test server)

### 🧪 2. UI Test Panel
**Interactive testing within the app interface**

```bash
# Start dev server first
npm run dev

# Then use keyboard shortcut or click button
Cmd/Ctrl+T  # Open test panel
```

**Features:**
- 🎯 8 predefined test objects with expected results
- 🚀 "Run All Tests" for comprehensive testing
- 🎨 Color-coded categories (food, beverage, chemical, etc.)
- 💭 Hover to see expected molecules
- 📊 Real-time console output with detailed results

### 💻 3. Command Line Tester
**Standalone script for detailed technical feedback**

```bash
# Start dev server in one terminal
npm run dev

# Run tests in another terminal
./test-molecular-ui.js              # Run all predefined tests
./test-molecular-ui.js -i           # Interactive mode - type any object
./test-molecular-ui.js --help       # Show help
```

**Output includes:**
- Detailed step-by-step pipeline execution
- Molecule counts and names
- SDF file generation status
- File accessibility verification
- Success/failure rates

## Chemical Test Data (Minimal Subset Approach)

Tests validate that expected chemicals are found as a **minimal subset** within the AI response. Additional compounds in the response are acceptable and expected.

### Basics (High Confidence Minimal Subset)
- **Water**: `O` (expect: water)
- **Ethanol**: `CCO` (expect: ethanol)  
- **Sodium Chloride**: `[Na+].[Cl-]` (expect: sodium, chloride)

### Beverages (Essential Components Only)
- **Red Wine**: (expect: ethanol, water) + may contain tartaric acid, glucose, tannins, etc.
- **Black Coffee**: (expect: water, caffeine) + may contain chlorogenic acid, etc.

### Biological (Most Confident Components)
- **Fresh Apple**: (expect: water) + may contain fructose, glucose, cellulose, malic acid, etc.

## Usage Recommendations

### Development Workflow
1. **Quick UI Testing**: Use test panel (`Cmd/Ctrl+T`) for rapid feedback
2. **Detailed Analysis**: Use command line tester for technical details
3. **CI/CD Integration**: Use `npm run test:pipeline` for automated testing

### Debugging Issues
1. **No molecules found**: Check AI analysis step in command line tester
2. **No visualization**: Verify SDF generation and file serving
3. **Poor quality**: Review SMILES validation in full pipeline test

### Performance Validation
- Full pipeline should complete within 30 seconds
- SDF generation success rate should be >70%
- File accessibility should be >70%

## Test Coverage

### What Gets Tested
✅ **Complete Pipeline**: Analysis → SDF → Visualization  
✅ **Quality Validation**: SMILES syntax, chemical accuracy  
✅ **Error Handling**: Malformed JSON, invalid SMILES  
✅ **Performance**: Response times, success rates  
✅ **File Management**: SDF creation, HTTP serving  

### What's NOT Tested
❌ 3DMol.js rendering (requires browser)  
❌ UI interactions (requires E2E tests)  
❌ Visual appearance (requires screenshot comparison)  

## Troubleshooting

### Common Issues

**"Skipping tests - no OPENAI_API_KEY"**
```bash
export OPENAI_API_KEY="your-key-here"
npm run test:pipeline
```

**"Failed to parse AI response"**
- This is expected for complex molecules
- Fallback handlers should provide backup data
- Check console for repair attempts

**"SDF generation failed"**
- Some SMILES may be invalid (expected)
- Check error count vs success count
- Ensure Python/RDKit environment is working

**"File not accessible"**
- Verify dev server is running on correct port
- Check `/sdf_files/` static file serving
- Ensure SDF directory exists and has files

## Example Output

```
   🧪 Testing full pipeline for: red wine
   Category: beverages
   Expected components: 2 (minimal subset)
   📊 Step 1: AI Analysis...
   ✅ Analysis complete: 6 molecules found
      1. Ethanol (CCO)
      2. Water (O)
      3. Tartaric acid (OC(C(O)C(O)=O)C(O)=O)
      4. Glucose (C(C(C(C(C(C=O)O)O)O)O)O)
      5. Resveratrol (C1=CC(=CC=C1C=CC2=CC(=CC(=C2)O)O)O)
      6. Malic acid (C(C(=O)O)C(C(=O)O)O)
   📁 Step 2: SDF Generation...
   ✅ SDF generation: 6 files created
   🔗 Step 3: File Accessibility...
   ✅ File access: 6/6 files accessible
   🎯 Step 4: Visualization Data...
   ✅ Visualization data ready for 6 molecules
   📈 Step 5: Quality Assessment...
   📋 Expected subset validation: 2/2 required chemicals found
      Expected subset: ethanol, water
      Found in response: ethanol, water
      Total response: 6 chemicals (may include additional valid compounds)
   🧬 SMILES quality: 100.0% (6/6)
   🎉 Pipeline completed successfully!
```

This comprehensive testing ensures your molecular visualization works correctly across all scenarios without requiring manual input each time!