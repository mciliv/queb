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

**Test Categories:**
- **Basics**: water, ethanol, sodium chloride (95% accuracy expected)
- **Beverages**: red wine, black coffee (75-80% accuracy expected)  
- **Biological**: fresh apple (65% accuracy expected)

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

## Chemical Test Data

All testing uses the same validated chemical compositions:

### Basics (High Accuracy Expected)
- **Water**: `O`
- **Ethanol**: `CCO`  
- **Sodium Chloride**: `[Na+].[Cl-]`

### Beverages (Realistic Composition)
- **Red Wine**: Ethanol, water, tartaric acid
- **Black Coffee**: Water, caffeine

### Biological (Complex but Realistic)
- **Fresh Apple**: Water, fructose, glucose, cellulose

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
   Expected components: 3
   📊 Step 1: AI Analysis...
   ✅ Analysis complete: 4 molecules found
      1. Ethanol (CCO)
      2. Water (O)
      3. Tartaric acid (OC(C(O)C(O)=O)C(O)=O)
      4. Glucose (C(C(C(C(C(C=O)O)O)O)O)O)
   📁 Step 2: SDF Generation...
   ✅ SDF generation: 4 files created
   🔗 Step 3: File Accessibility...
   ✅ File access: 4/4 files accessible
   🎯 Step 4: Visualization Data...
   ✅ Visualization data ready for 4 molecules
   📈 Step 5: Quality Assessment...
   📋 Required chemicals found: 3/3
   🧬 SMILES quality: 100.0% (4/4)
   🎉 Pipeline completed successfully!
```

This comprehensive testing ensures your molecular visualization works correctly across all scenarios without requiring manual input each time!