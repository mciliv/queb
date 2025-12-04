/**
 * Molecular Analysis Schemas
 */

const { z } = require('zod');

// SMILES array validation
const smilesArray = z.object({
  smiles: z.array(z.string()).min(1, 'At least one SMILES string is required')
});

// Object molecules request validation
const objectMoleculesRequest = z.object({
  object: z.string().min(1, 'Object description is required')
});

// SDF generation request validation
const sdfGenerationRequest = z.object({
  smiles: z.array(z.string()).min(1, 'At least one SMILES string is required'),
  overwrite: z.boolean().optional().default(false)
});

// Image molecule analysis request
const ImageMoleculeSchema = z.object({
  image: z.string().min(1, 'Image data is required'),
  format: z.enum(['base64', 'url']).optional().default('base64')
});

// Text molecule analysis request
const TextMoleculeSchema = z.object({
  text: z.string().min(1, 'Text input is required'),
  minMoleculeCount: z.number().min(1).optional().default(1)
});

// SDF generation schema
const SdfGenerationSchema = z.object({
  smiles: z.array(z.string()).min(1, 'At least one SMILES string is required'),
  overwrite: z.boolean().optional().default(false)
});

module.exports = {
  smilesArray,
  objectMoleculesRequest,
  sdfGenerationRequest,
  ImageMoleculeSchema,
  TextMoleculeSchema,
  SdfGenerationSchema
};
