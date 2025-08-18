// Google Cloud Functions entry point for Mol
// =========================================

const express = require("express");
const cors = require("cors");
const fs = require("fs");
const path = require("path");
const Structuralizer = require("../services/StructurePredictor");
const MolecularProcessor = require("../services/molecular-processor");
const {
  ImageMoleculeSchema,
  TextMoleculeSchema,
  SdfGenerationSchema,
} = require("../schemas/schemas");

// Initialize Express app
const app = express();

// Initialize modules
const structuralizer = new Structuralizer(process.env.OPENAI_API_KEY);
const molecularProcessor = new MolecularProcessor();

// ==================== MIDDLEWARE ====================
app.use(cors());
app.use(express.json({ limit: "50mb" }));

// ==================== ROUTES ====================

// Shared text analysis handler (single implementation, multiple aliases)
const handleTextAnalysis = async (req, res) => {
  try {
    const validation = TextMoleculeSchema.safeParse(req.body);
    if (!validation.success) {
      return res.status(400).json({
        error: "Invalid input data",
        details: validation.error.issues,
      });
    }
    const { object } = req.body;
    if (!object || typeof object !== 'string' || object.trim().length === 0) {
      return res.status(400).json({ error: "No object description provided" });
    }
    const result = await structuralizer.analyzeText(object.trim());
    res.json({ output: result });
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
};

// TEXT FIRST (primary + aliases)
app.post("/structuralize-text", handleTextAnalysis);
app.post("/analyze-text", handleTextAnalysis);
app.post("/structures-from-text", handleTextAnalysis);
app.post("/object-molecules", handleTextAnalysis);

// Image analysis route
app.post("/image-molecules", async (req, res) => {
  try {
    const validation = ImageMoleculeSchema.safeParse(req.body);
    if (!validation.success) {
      return res.status(400).json({
        error: "Invalid input data",
        details: validation.error.issues,
      });
    }

    const {
      imageBase64,
      croppedImageBase64,
      x,
      y,
      cropMiddleX,
      cropMiddleY,
      cropSize,
    } = req.body;

    if (!imageBase64) {
      return res.status(400).json({ error: "No image data provided" });
    }

    const result = await structuralizer.analyzeImage(
      imageBase64,
      croppedImageBase64,
      x,
      y,
      cropMiddleX,
      cropMiddleY,
      cropSize,
    );
    res.json({ output: result });
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
});

// SDF generation route
app.post("/generate-sdfs", async (req, res) => {
  try {
    // Validate input schema
    const validation = SdfGenerationSchema.safeParse(req.body);
    if (!validation.success) {
      return res.status(400).json({
        error: "Invalid input data",
        details: validation.error.issues,
      });
    }

    const { smiles, overwrite = false } = req.body;

    if (!smiles || !Array.isArray(smiles)) {
      return res.status(400).json({ error: "smiles array is required" });
    }

    const result = await molecularProcessor.processSmiles(smiles, overwrite);

    res.json({
      sdfPaths: result.sdfPaths,
      errors: result.errors,
      skipped: result.skipped,
      message: `Generated ${result.sdfPaths.length} 3D structures from ${smiles.length} SMILES`,
    });
  } catch (error) {

    res.status(500).json({ error: error.message });
  }
});

// Health check route
app.get("/", (req, res) => {
  res.json({
    status: "ok",
    message: "Mol Molecular Analysis API",
    version: "1.0.0",
  });
});

// Export for Google Cloud Functions
exports.main = app;
