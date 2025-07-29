/**
 * MOLECULAR PROCESSOR SERVICE
 * Purpose: Convert SMILES notation to 3D molecular structures (SDF files)
 * Features: SMILES validation, SDF generation, file management, error handling
 * 
 * Core Responsibilities:
 * - Validate SMILES notation format
 * - Generate 3D molecular structures from SMILES
 * - Manage SDF file storage and retrieval
 * - Handle batch processing of multiple molecules
 * - Interface with Python chemistry processing scripts
 */

const { spawn } = require("child_process");
const fs = require("fs");
const path = require("path");
const { FILE_CONFIG, CHEMISTRY_CONFIG } = require("../config/constants");

class MolecularProcessor {
  constructor(sdfDir = FILE_CONFIG.SDF_DIRECTORY) {
    this.sdfDir = path.join(__dirname, "..", "..", sdfDir);
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
        // Basic SMILES validation
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

  // Basic SMILES format validation using centralized patterns
  isValidSmiles(smiles) {
    if (!smiles || typeof smiles !== 'string') return false;
    if (CHEMISTRY_CONFIG.INVALID_SMILES.includes(smiles.trim())) return false;
    
    // Skip obvious molecular formulas using centralized pattern
    if (CHEMISTRY_CONFIG.MOLECULAR_FORMULA_PATTERN.test(smiles.replace(/\s/g, ''))) {
      return false; // Likely molecular formula like H2O, CaCO3, etc.
    }
    
    return true;
  }

  async generateSDF(smiles, overwrite = false) {
    // Check if SDF already exists
    if (!overwrite) {
      const existingPath = this.findExistingSdfFile(smiles);
      if (existingPath) {
        console.log(`Using existing file: ${smiles} → ${existingPath}`);
        return existingPath;
      }
    }

    // Generate from SMILES
    try {
      console.log(`Generating SMILES structure for: ${smiles}`);
      const sdfPath = await this.generateSmilesSDF(smiles);
      if (sdfPath) return sdfPath;
    } catch (error) {
      console.error(`SMILES generation failed for ${smiles}:`, error.message);
      console.error(`Full error:`, error);
      throw error;
    }

    return null;
  }

  async generateSmilesSDF(chemical) {
    return new Promise((resolve, reject) => {
      const pythonProcess = spawn("python", [
        path.join(__dirname, "..", "..", "chemistry", "processors", "sdf.py"),
        chemical,
        "--dir",
        this.sdfDir,
      ], {
        cwd: path.join(__dirname, "..", "..")
      });

      let output = "";
      pythonProcess.stdout.on("data", (data) => {
        output += data.toString();
      });

      pythonProcess.stderr.on("data", (data) => {
        console.error(`Python Error: ${data.toString()}`);
      });

      pythonProcess.on("close", (code) => {
        if (code === 0) {
          const sdfPath = this.findExistingSdfFile(chemical);
          if (sdfPath) {
            console.log(
              `Successfully generated SMILES structure: ${chemical} → ${sdfPath}`,
            );
            resolve(sdfPath);
          } else {
            reject(
              new Error(
                `SMILES succeeded but couldn't find SDF file for ${chemical}`,
              ),
            );
          }
        } else {
          console.error(`Python process failed with exit code: ${code}`);
          console.error(`Working directory: ${process.cwd()}`);
          console.error(`Python command: python sdf.py ${chemical} --dir ${this.sdfDir}`);
          console.error(`SDF directory exists: ${fs.existsSync(this.sdfDir)}`);
          reject(new Error(`SMILES generation failed for ${chemical} (exit code: ${code})`));
        }
      });
    });
  }

  findExistingSdfFile(smiles) {
    const possibleFilenames = [
      `${smiles}.sdf`,
      `${smiles.replace(/[^a-zA-Z0-9]/g, "_")}.sdf`,
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
