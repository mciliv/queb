# Chemistry Processing

This directory contains chemical data processing utilities for the molecular space analyzer web application.

## Components

### `/processors/`
- **`sdf.py`** - SMILES to SDF (Structure Data File) conversion
  - Converts SMILES strings to 3D molecular structures
  - Used by the web app for 3D visualization with 3DMol.js
  - Generates files for the molecular viewer components

## Purpose

Provides the chemical format conversion needed for the web application's molecular visualization features. The SDF files generated here are consumed by the frontend's 3D molecular viewer to display interactive molecular structures.

## Dependencies

- RDKit - Chemical informatics toolkit
- Logging - For conversion process tracking

## Usage

This module is called by the backend molecular processor service when generating 3D structures from SMILES data returned by the AI analysis.