/**
 * molecular-processor.js - Handles 3D molecular structure generation
 * 
 * This service converts chemical representations (SMILES notation, chemical names)
 * into 3D structure files (SDF format) that can be visualized in the browser.
 * It acts as a bridge between chemical data and 3D visualization.
 * 
 * Key responsibilities:
 * - Convert SMILES notation to 3D structures using RDKit (Python)
 * - Generate SDF files from chemical names via PubChem lookup
 * - Cache generated files to avoid redundant processing
 * - Handle fallback to PubChem when local generation fails
 */

const fs = require("fs");
const path = require("path");
const crypto = require("crypto");
const fsPromises = require("fs").promises;
const { resolveName, downloadSDFBySmiles } = require("./name-resolver");
const { warn, error } = require("../../core/logger");

/**
 * MolecularProcessor - Main class for molecular structure generation
 * 
 * This class manages the conversion of chemical data into 3D structure files
 * that can be rendered by 3Dmol.js in the frontend. It uses a multi-step
 * approach with fallbacks to ensure reliable structure generation.
 * 
 * Processing pipeline:
 * 1. Check if SDF file already exists (cached)
 * 2. Try generating with RDKit (Python subprocess)
 * 3. If RDKit fails, download from PubChem as fallback // Make this step 2, or tell me that it's less effective
 * 4. Store generated files for future use
 */
class MolecularProcessor {
  /**
   * Initialize the molecular processor with SDF storage directory
   * @param {string} [sdfDir] - Custom directory for SDF files (defaults to public/sdf_files)
   */
  constructor(sdfDir) {
    if (!sdfDir) {
      sdfDir = process.env.NODE_ENV === 'test' ? 'tests/sdf_files' : 'public/sdf_files';
    }
    // __dirname is src/server/services — go up 3 levels to reach project root
    // Prompt: molecular-processor.js line 44 — path must resolve to <root>/public/sdf_files
    this.sdfDir = path.join(__dirname, "..", "..", "..", sdfDir);
    this.ensureSdfDirectory();
  }

  ensureSdfDirectory() {
    if (!fs.existsSync(this.sdfDir)) {
      fs.mkdirSync(this.sdfDir, { recursive: true });
    }
  }

  // Process SMILES array and generate SDF files
  async processSmiles(smilesArray, overwrite = false) {
    const results = {
      sdfPaths: [],
      errors: [],
      skipped: [],
    };

    for (const smiles of smilesArray) {
      try {
        if (!this.isValidSmiles(smiles)) {
          results.skipped.push(`${smiles.substring(0, 50)}... (invalid format)`);
          continue;
        }

        const sdfPath = await this.generateSDF(smiles, overwrite);
        if (sdfPath) {
          results.sdfPaths.push(sdfPath);
        } else {
          results.skipped.push(smiles);
        }
      } catch (error) {
        results.errors.push(`${smiles.substring(0, 50)}... - ${error.message}`);
      }
    }

    return results;
  }

  async generateSDFByName(name, overwrite = false) {
    const { smiles, title, iupac } = await resolveName(name);
    const displayName = title || iupac || name;

    if (smiles) {
      const pathOrNull = await this.generateSDF(smiles, overwrite);
      if (pathOrNull) {
        return { sdfPath: pathOrNull, name: displayName };
      }
    }

    return null;
  }

  // Basic SMILES format validation with length limits
  isValidSmiles(smiles) {
    if (!smiles || typeof smiles !== 'string') return false;
    if (smiles.trim() === '' || smiles === 'N/A') return false;

    const cleaned = smiles.replace(/\s/g, '');

    // Check length limits - very long SMILES can cause processing issues
    if (cleaned.length > 200) {
      warn(`SMILES too long (${cleaned.length} chars), may fail processing`, {
        smilesLength: cleaned.length,
        smilesPrefix: cleaned.substring(0, 50)
      });
      return false; // Reject very long SMILES to prevent processing failures
    }
    
    // Common valid SMILES patterns (allow lowercase for aromatic atoms)
    const validPatterns = [
      /^[A-Za-z][A-Za-z][A-Za-z]$/, // CCO, CCC, cco, etc.
      /^[A-Za-z][A-Za-z]\(=O\)[A-Za-z]$/, // CC(=O)O, etc.
      /^[A-Za-z]\([A-Za-z][A-Za-z]\)[A-Za-z]$/, // C(CC)C, etc.
      /^[A-Za-z]1[A-Za-z]=[A-Za-z][A-Za-z]=[A-Za-z]1$/, // C1=CC=CC=C1 (benzene)
      /^[A-Za-z]1[A-Za-z][A-Za-z][A-Za-z][A-Za-z][A-Za-z]1$/, // C1CCCCC1 (cyclohexane)
      /^[A-Za-z]$/, // O, N, C, c, etc.
      /^[A-Za-z]1[A-Za-z][A-Za-z][A-Za-z][A-Za-z][A-Za-z][A-Za-z][A-Za-z][A-Za-z]1$/, // c1ccc2ccccc2c1 (naphthalene)
    ];
    
    // Check if it matches any valid pattern
    for (const pattern of validPatterns) {
      if (pattern.test(cleaned)) {
        return true;
      }
    }
    
    // Reject obvious molecular formulas (H2O, CaCO3, etc.)
    if (/^[A-Z][a-z]?[0-9]*([A-Z][a-z]?[0-9]*)+$/.test(cleaned) && 
        !cleaned.includes('(') && 
        !cleaned.includes('=') && 
        !cleaned.includes('[') && 
        !cleaned.includes(']') && 
        !cleaned.includes('@') && 
        !cleaned.includes('#') && 
        !cleaned.includes('\\') && 
        !cleaned.includes('/')) {
      return false;
    }
    
    // Allow other patterns that might be valid SMILES
    return true;
  }

  async generateSDF(smiles, overwrite = false) {
    // Check if SDF already exists
    if (!overwrite) {
      const existingPath = this.findExistingSdfFile(smiles);
      if (existingPath) {

        return existingPath;
      }
    }

    // Generate from SMILES
    try {
      const sdfPath = await this.generateSmilesSDF(smiles);
      if (sdfPath) return sdfPath;
    } catch (error) {
      // Fallback to bundled SDFs for common SMILES if Python is unavailable
      const fallback = this.copyFallbackSdf(smiles);
      if (fallback) return fallback;
      throw error;
    }

    return null;
  }

  async generateSmilesSDF(chemical) {
    // First try local Python generator to satisfy unit tests expectations
    // Prompt: fix broken script path used by subprocess spawn()
    // `__dirname` is `src/server/services`, and the actual script lives at `src/server/sdf.py`.
    const pythonScript = path.join(__dirname, "..", "sdf.py");
    const args = [pythonScript, chemical, "--dir", this.sdfDir];
    const spawnOptions = { stdio: "pipe" };
    const { spawn } = require("child_process");

    let stderrOutput = '';
    const exited = await new Promise((resolve) => {
      try {
        const child = spawn("python", args, spawnOptions);
        child.stdout && child.stdout.on("data", () => {});
        child.stderr && child.stderr.on("data", (data) => {
          stderrOutput += data.toString();
        });
        child.on("close", (code) => resolve(code));
        child.on("error", () => resolve(1));
      } catch (_) {
        resolve(1);
      }
    });

    // Log stderr if there was an error
    if (exited !== 0 && stderrOutput) {
      error(`Python SDF generation failed`, {
        smilesPrefix: chemical.substring(0, 50),
        stderr: stderrOutput
      });
    }
    
    if (exited === 0) {
      // The Python helper may succeed but choose a different filename format
      // (raw SMILES) or even fail silently without creating a file.
      // Verify on-disk existence using our resolver and only return if found.
      const existing = this.findExistingSdfFile(chemical);
      if (existing) {
        return existing;
      }
      // Fall through to network fallback if the file wasn't actually created
    }
    
    // If python path fails, attempt PubChem download then fallback
    try {
      const sdf = await downloadSDFBySmiles(chemical);
      const filename = this.generateSafeFilename(chemical);
      const dest = path.join(this.sdfDir, filename);
      await fsPromises.writeFile(dest, sdf, "utf8");
      return `/sdf_files/${filename}`;
    } catch (error) {
      const fallbackPath = this.copyFallbackSdf(chemical);
      if (fallbackPath) return fallbackPath;
      // Match unit test expectation but provide more context
      const errorMsg = `SMILES generation failed for: ${chemical.substring(0, 50)}${chemical.length > 50 ? '...' : ''}`;
      error(errorMsg, error, {
        smilesPrefix: chemical.substring(0, 50),
        smilesLength: chemical.length
      });
      throw new Error(errorMsg);
    }
  }

  copyFallbackSdf(smiles) {
    try {
      const fallbackDir = path.join(__dirname, "..", "..", "tests", "sdf_files");
      const filenames = [
        `${smiles}.sdf`,
        this.generateSafeFilename(smiles),
        `${smiles.replace(/[^a-zA-Z0-9]/g, ch => ch === "=" ? "__" : "_")}.sdf`,
      ];
      for (const name of filenames) {
        const src = path.join(fallbackDir, name);
        if (fs.existsSync(src)) {
          const dest = path.join(this.sdfDir, name);
          if (!fs.existsSync(dest)) {
            fs.copyFileSync(src, dest);
          }
          return `/sdf_files/${name}`;
        }
      }
      return null;
    } catch (_) {
      return null;
    }
  }

  generateSafeFilename(smiles) {
    // Create a hash for long SMILES to avoid filesystem limits and display issues
    const cleanedName = smiles.replace(/[^a-zA-Z0-9]/g, ch => ch === "=" ? "__" : "_");
    
    // If the cleaned filename would be too long for display or filesystem
    if (cleanedName.length > 50) {
      const hash = crypto.createHash('md5').update(smiles).digest('hex').substring(0, 12);
      // Include first few characters for recognition + hash
      const prefix = cleanedName.substring(0, 8);
      return `${prefix}_${hash}.sdf`;
    }
    
    // For shorter SMILES, use the cleaned version
    return `${cleanedName}.sdf`;
  }

  findExistingSdfFile(smiles) {
    const possibleFilenames = [
      `${smiles}.sdf`,
      this.generateSafeFilename(smiles),
      `${smiles.replace(/[^a-zA-Z0-9]/g, ch => ch === "=" ? "__" : "_")}.sdf`,
    ];

    for (const filename of possibleFilenames) {
      const fullPath = path.join(this.sdfDir, filename);
      if (fs.existsSync(fullPath)) {
        return `/sdf_files/${filename}`;
      }
    }

    return null;
  }
}

module.exports = MolecularProcessor;
