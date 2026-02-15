const { z } = require("zod");

// Chemical representations - focus on SMILES constituents
const CHEMICAL_REPRESENTATIONS = {
  smiles:
    "Array of SMILES notation for constituent molecules (e.g., ['CCO', 'O'] for alcoholic beverage, ['C1=CC=CC=C1'] for benzene)",
};

// Schema for image molecule analysis
const ImageMoleculeSchema = z.object({
  imageBase64: z.string().min(1, "Image data is required"),
  croppedImageBase64: z.string().optional(),
  x: z.number().min(0).max(1000).optional(),
  y: z.number().min(0).max(1000).optional(),
  cropMiddleX: z.number().min(0).max(1000).optional(),
  cropMiddleY: z.number().min(0).max(1000).optional(),
  cropSize: z.number().min(10).max(500).optional(),
});

// Schema for text molecule analysis
const TextMoleculeSchema = z.object({
  object: z.string()
    .min(1, "Object description is required")
    .refine(
      (value) => {
        // Basic validation for physical objects
        const text = value.toLowerCase().trim();
        
        // Reject obviously non-physical things
        const nonPhysicalPatterns = [
          /^(love|hate|happy|sad|angry|joy|fear|hope|dream|idea|thought|feeling|emotion)/,
          /^(running|walking|talking|thinking|sleeping|eating|drinking)$/,
          /^[a-z]{1,2}$/,  // Single letters or very short nonsense
          /^[^a-z]*$/,     // Only numbers/symbols
          /^(asdf|qwerty|test|random|nothing|something|anything|everything)$/i
        ];
        
        return !nonPhysicalPatterns.some(pattern => pattern.test(text));
      },
      {
        message: "Please describe a real, physical object (food, materials, plants, etc.)"
      }
    )
});

// Schema for SDF generation
const SdfGenerationSchema = z.object({
  smiles: z.array(z.string()).min(1, "At least one SMILES string is required"),
  overwrite: z.boolean().optional(),
});

// Schema for object identification response
const ObjectIdentificationSchema = z.object({
  object: z.string(),
  chemicals: z.array(
    z.object({
      name: z.string(),
      smiles: z.string(),
    }),
  ),
});

module.exports = {
  ImageMoleculeSchema,
  TextMoleculeSchema,
  SdfGenerationSchema,
  ObjectIdentificationSchema,
};
