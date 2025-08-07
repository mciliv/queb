# Molecular Docking Research Archive

This directory contains the molecular docking and computational chemistry research code from earlier versions of this project (commits before molecular space analyzer web app).

## Archived Components

### `/dock/` - Molecular Docking Pipeline
- **`chem.py`** - Chemical structure manipulation and PDB processing
- **`dock.py`** - Main docking workflow with GNINA integration  
- **`rcsb_pdb.py`** - RCSB Protein Data Bank retrieval and processing

### Conversion Processors
- **`crystal.py`** - Crystal structure generation for minerals
- **`uniprot.py`** - UniProt protein sequence processing
- **`usdz.py`** - SMILES to USDZ 3D format conversion

### `/models/` - Machine Learning Models
- Research models for molecular property prediction

## Separation Rationale

This code was separated from the main molecular space analyzer web application because:

1. **Different Purpose**: Research/academic molecular docking vs. web-based molecular identification
2. **Different Dependencies**: Heavy computational chemistry tools vs. lightweight web stack
3. **Memory Focus**: Cleaner context for current web application development

## Original Functionality

The archived code provided:
- Automated molecular docking with GNINA
- PDB structure manipulation and analysis
- Ligand-protein interaction studies
- Multiple output format conversions (SDF, USDZ, etc.)
- Integration with RCSB PDB and UniProt databases

## Current Usage

The main project now only uses `molecular-conversion/processors/sdf.py` for SMILES → SDF conversion in the web application.

## Historical Context

This separation was made on [current date] to focus the main repository on the molecular space analyzer web application while preserving the valuable research code for potential future use.