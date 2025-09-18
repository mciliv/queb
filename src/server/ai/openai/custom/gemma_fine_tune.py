#!/usr/bin/env python3
"""
Fine-tune Gemma for molecular analysis tasks.
Supports SMILES notation, protein analysis, and molecular property prediction.
"""

import os
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field
import torch
from torch.utils.data import Dataset, DataLoader
from transformers import (
    AutoTokenizer, 
    AutoModelForCausalLM,
    TrainingArguments,
    Trainer,
    DataCollatorForLanguageModeling
)
from peft import LoraConfig, get_peft_model, TaskType
import pandas as pd
from datasets import Dataset as HFDataset

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class MolecularTrainingConfig:
    """Configuration for molecular fine-tuning."""
    model_name: str = "google/gemma-2b"
    output_dir: str = "./gemma-mol-finetuned"
    data_dir: str = "./backend/molecular-docking-research"
    max_length: int = 512
    batch_size: int = 4
    learning_rate: float = 2e-4
    num_epochs: int = 3
    warmup_steps: int = 100
    save_steps: int = 500
    logging_steps: int = 100
    
    # LoRA configuration
    lora_r: int = 16
    lora_alpha: int = 32
    lora_dropout: float = 0.1
    
    # Training data types
    include_smiles: bool = True
    include_protein_analysis: bool = True
    include_docking_results: bool = True

class MolecularDataset(Dataset):
    """Dataset for molecular fine-tuning with SMILES and protein data."""
    
    def __init__(self, config: MolecularTrainingConfig, tokenizer):
        self.config = config
        self.tokenizer = tokenizer
        self.data = self._load_molecular_data()
        
    def _load_molecular_data(self) -> List[Dict[str, str]]:
        """Load and format molecular data for training."""
        # Try to load extracted molecular data first
        extracted_data_file = "molecular_training_data.json"
        if Path(extracted_data_file).exists():
            logger.info(f"Loading extracted molecular data from {extracted_data_file}")
            with open(extracted_data_file, 'r') as f:
                return json.load(f)
        
        # Fallback to synthetic data
        logger.info("Using synthetic molecular data")
        data = []
        
        if self.config.include_smiles:
            data.extend(self._load_smiles_data())
            
        if self.config.include_protein_analysis:
            data.extend(self._load_protein_data())
            
        if self.config.include_docking_results:
            data.extend(self._load_docking_data())
            
        return data
    
    def _load_smiles_data(self) -> List[Dict[str, str]]:
        """Generate training data from SMILES notation."""
        smiles_data = []
        
        # Training examples for SMILES interpretation
        examples = [
            {
                "instruction": "Convert this molecule name to SMILES notation:",
                "input": "caffeine",
                "output": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
            },
            {
                "instruction": "Describe the molecular properties of this SMILES:",
                "input": "CCO",
                "output": "This is ethanol (C2H6O), a simple alcohol with a hydroxyl group (-OH) attached to a two-carbon chain. It has hydrogen bonding capability and is polar."
            },
            {
                "instruction": "Identify functional groups in this SMILES:",
                "input": "CC(=O)O",
                "output": "This molecule contains a carboxyl group (-COOH), making it acetic acid. The functional group includes both a carbonyl (C=O) and hydroxyl (OH) group."
            }
        ]
        
        for example in examples:
            prompt = f"Instruction: {example['instruction']}\nInput: {example['input']}\nOutput: {example['output']}"
            smiles_data.append({"text": prompt})
            
        return smiles_data
    
    def _load_protein_data(self) -> List[Dict[str, str]]:
        """Generate training data from protein analysis."""
        protein_data = []
        
        examples = [
            {
                "instruction": "Analyze this PDB structure:",
                "input": "ATOM      1  N   ALA A   1      20.154  16.967  23.986  1.00 20.00           N",
                "output": "This is an alanine residue (ALA) in chain A, position 1. The nitrogen atom coordinates are (20.154, 16.967, 23.986) with occupancy 1.00 and B-factor 20.00."
            },
            {
                "instruction": "Explain protein-ligand interaction:",
                "input": "docking score: -8.5 kcal/mol",
                "output": "A docking score of -8.5 kcal/mol indicates strong binding affinity between the protein and ligand. This suggests favorable interactions and good complementarity."
            }
        ]
        
        for example in examples:
            prompt = f"Instruction: {example['instruction']}\nInput: {example['input']}\nOutput: {example['output']}"
            protein_data.append({"text": prompt})
            
        return protein_data
    
    def _load_docking_data(self) -> List[Dict[str, str]]:
        """Generate training data from molecular docking results."""
        docking_data = []
        
        examples = [
            {
                "instruction": "Interpret molecular docking results:",
                "input": "Receptor: 2nnq_protein.pdb, Ligand: caffeine.sdf, Score: -7.2",
                "output": "The caffeine ligand shows good binding to the 2nnq protein receptor with a score of -7.2 kcal/mol, indicating favorable molecular interactions."
            }
        ]
        
        for example in examples:
            prompt = f"Instruction: {example['instruction']}\nInput: {example['input']}\nOutput: {example['output']}"
            docking_data.append({"text": prompt})
            
        return docking_data
    
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        item = self.data[idx]
        encoding = self.tokenizer(
            item["text"],
            truncation=True,
            padding="max_length",
            max_length=self.config.max_length,
            return_tensors="pt"
        )
        
        return {
            "input_ids": encoding["input_ids"].flatten(),
            "attention_mask": encoding["attention_mask"].flatten(),
            "labels": encoding["input_ids"].flatten()
        }

def setup_model_and_tokenizer(config: MolecularTrainingConfig):
    """Initialize Gemma model and tokenizer with LoRA."""
    # Load tokenizer
    tokenizer = AutoTokenizer.from_pretrained(config.model_name)
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token
    
    # Load model
    model = AutoModelForCausalLM.from_pretrained(
        config.model_name,
        torch_dtype=torch.float16 if torch.cuda.is_available() else torch.float32,
        device_map="auto" if torch.cuda.is_available() else None,
        trust_remote_code=True
    )
    
    # Configure LoRA
    lora_config = LoraConfig(
        task_type=TaskType.CAUSAL_LM,
        r=config.lora_r,
        lora_alpha=config.lora_alpha,
        lora_dropout=config.lora_dropout,
        target_modules=["q_proj", "k_proj", "v_proj", "o_proj"]
    )
    
    model = get_peft_model(model, lora_config)
    model.print_trainable_parameters()
    
    return model, tokenizer

def train_gemma_molecular(config: MolecularTrainingConfig):
    """Fine-tune Gemma for molecular analysis tasks."""
    logger.info("Starting Gemma molecular fine-tuning")
    
    # Setup model and tokenizer
    model, tokenizer = setup_model_and_tokenizer(config)
    
    # Create dataset
    dataset = MolecularDataset(config, tokenizer)
    logger.info(f"Created dataset with {len(dataset)} examples")
    
    # Convert to HuggingFace dataset
    hf_dataset = HFDataset.from_list(dataset.data)
    
    def tokenize_function(examples):
        return tokenizer(
            examples["text"],
            truncation=True,
            padding="max_length",
            max_length=config.max_length
        )
    
    tokenized_dataset = hf_dataset.map(tokenize_function, batched=True)
    
    # Training arguments
    training_args = TrainingArguments(
        output_dir=config.output_dir,
        num_train_epochs=config.num_epochs,
        per_device_train_batch_size=config.batch_size,
        gradient_accumulation_steps=2,
        warmup_steps=config.warmup_steps,
        learning_rate=config.learning_rate,
        fp16=torch.cuda.is_available(),
        logging_steps=config.logging_steps,
        save_steps=config.save_steps,
        evaluation_strategy="no",
        save_strategy="steps",
        load_best_model_at_end=False,
        report_to=None,
        remove_unused_columns=False,
    )
    
    # Data collator
    data_collator = DataCollatorForLanguageModeling(
        tokenizer=tokenizer,
        mlm=False,
        pad_to_multiple_of=8 if torch.cuda.is_available() else None,
    )
    
    # Trainer
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=tokenized_dataset,
        data_collator=data_collator,
        tokenizer=tokenizer,
    )
    
    # Train
    logger.info("Starting training...")
    trainer.train()
    
    # Save model
    trainer.save_model()
    tokenizer.save_pretrained(config.output_dir)
    
    logger.info(f"Training completed. Model saved to {config.output_dir}")

def test_molecular_inference(model_path: str, test_prompts: List[str]):
    """Test the fine-tuned model on molecular tasks."""
    logger.info("Testing molecular inference")
    
    tokenizer = AutoTokenizer.from_pretrained(model_path)
    model = AutoModelForCausalLM.from_pretrained(
        model_path,
        torch_dtype=torch.float16 if torch.cuda.is_available() else torch.float32,
        device_map="auto" if torch.cuda.is_available() else None,
    )
    
    for prompt in test_prompts:
        inputs = tokenizer(prompt, return_tensors="pt")
        
        with torch.no_grad():
            outputs = model.generate(
                inputs.input_ids,
                max_new_tokens=100,
                temperature=0.7,
                do_sample=True,
                pad_token_id=tokenizer.eos_token_id
            )
        
        response = tokenizer.decode(outputs[0], skip_special_tokens=True)
        logger.info(f"Prompt: {prompt}")
        logger.info(f"Response: {response[len(prompt):]}")
        print("-" * 50)

def main():
    """Main training function."""
    config = MolecularTrainingConfig()
    
    # Create output directory
    os.makedirs(config.output_dir, exist_ok=True)
    
    # Save configuration
    with open(Path(config.output_dir) / "training_config.json", "w") as f:
        json.dump(config.__dict__, f, indent=2)
    
    # Train model
    train_gemma_molecular(config)
    
    # Test inference
    test_prompts = [
        "Instruction: Convert this molecule name to SMILES notation:\nInput: aspirin\nOutput:",
        "Instruction: Analyze this molecular docking score:\nInput: -9.2 kcal/mol\nOutput:",
        "Instruction: Describe this SMILES structure:\nInput: CC(C)(C)OC(=O)NC1=CC=C(C=C1)O\nOutput:"
    ]
    
    test_molecular_inference(config.output_dir, test_prompts)

if __name__ == "__main__":
    main()