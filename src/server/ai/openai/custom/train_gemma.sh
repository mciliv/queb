#!/bin/bash
set -e

echo "🧬 Starting Gemma Fine-tuning for Molecular Analysis"

# Check if Python 3.8+ is available
if ! python3 --version | grep -qE "3\.[89]|3\.1[0-9]"; then
    echo "❌ Python 3.8+ required for transformers"
    exit 1
fi

# Install dependencies if needed
if ! python3 -c "import torch" 2>/dev/null; then
    echo "📦 Installing PyTorch and dependencies..."
    pip3 install -r requirements.txt
fi

# Extract molecular data from codebase
echo "🔬 Extracting molecular data from existing workflows..."
python3 extract_molecular_data.py

# Check if GPU is available
if python3 -c "import torch; print(torch.cuda.is_available())" | grep -q "True"; then
    echo "🚀 GPU detected - using CUDA acceleration"
    export CUDA_VISIBLE_DEVICES=0
else
    echo "💻 Using CPU training (will be slower)"
fi

# Start fine-tuning
echo "🤖 Starting Gemma fine-tuning..."
python3 gemma_fine_tune.py

echo "✅ Gemma fine-tuning completed!"
echo "📁 Model saved to: ./gemma-mol-finetuned/"

# Test the fine-tuned model
echo "🧪 Testing fine-tuned model..."
python3 -c "
from gemma_fine_tune import test_molecular_inference

test_prompts = [
    'Instruction: Convert this molecule name to SMILES notation\nInput: water\nOutput:',
    'Instruction: Is this a valid SMILES structure?\nInput: CCO\nOutput:',
    'Instruction: Analyze the molecular structure of this SMILES\nInput: c1ccccc1\nOutput:'
]

test_molecular_inference('./gemma-mol-finetuned', test_prompts)
"

echo "🎉 Gemma fine-tuning for mol project complete!"