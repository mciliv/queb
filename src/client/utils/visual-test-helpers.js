/**
 * Visual test utilities - Pure functions for testable logic
 */

/**
 * Sanitizes SMILES notation for use in filenames
 * @param {string} smiles - SMILES notation string
 * @returns {string} Sanitized filename-safe string
 */
export function sanitizeSmiles(smiles) {
  return smiles.replace(/[^a-zA-Z0-9]/g, ch => ch === '=' ? '__' : '_');
}

/**
 * Normalizes paths by ensuring they start with a single forward slash
 * @param {string} path - File path to normalize
 * @returns {string} Normalized path
 */
export function normalizePath(path) {
  return path.replace(/^\/*/, '/');
}

/**
 * Resolves the best SDF path for a given SMILES string
 * Prefers server-returned paths, falls back to deterministic naming
 * 
 * @param {string} smiles - SMILES notation
 * @param {number} index - Index in the SMILES array
 * @param {string[]} returnedPaths - Array of paths returned from server
 * @returns {string} Resolved SDF file path
 */
export function resolveSdfPath(smiles, index, returnedPaths = []) {
  const deterministic = `/sdf_files/${sanitizeSmiles(smiles)}.sdf`;
  const byIndex = returnedPaths[index] || null;
  
  // Build a set of normalized returned paths for quick lookup
  const returnedFileSet = new Set(
    returnedPaths.map(p => normalizePath(p))
  );
  
  // Prefer deterministic path if server returned it, otherwise use index-based
  const normalizedDeterministic = normalizePath(deterministic);
  if (returnedFileSet.has(normalizedDeterministic)) {
    return deterministic;
  }
  
  return byIndex || deterministic;
}

/**
 * Creates a viewer object for molecular visualization
 * 
 * @param {string} smiles - SMILES notation
 * @param {string} sdfPath - Path to SDF file
 * @param {Object} nameMap - Map of SMILES to display names
 * @returns {Object} Viewer configuration object
 */
export function createViewerObject(smiles, sdfPath, nameMap = {}) {
  return {
    name: nameMap[smiles] || smiles,
    sdfData: sdfPath ? `file://${sdfPath}` : null,
    smiles: smiles
  };
}

/**
 * Processes a SMILES array into viewer objects
 * Main orchestration function for visual test setup
 * 
 * @param {string[]} smilesArray - Array of SMILES notations
 * @param {string[]} returnedPaths - Paths returned from SDF generation
 * @param {Object} nameMap - Map of SMILES to display names
 * @returns {Object[]} Array of viewer configuration objects
 */
export function processSmilesToViewers(smilesArray, returnedPaths = [], nameMap = {}) {
  return smilesArray.map((smiles, idx) => {
    const sdfPath = resolveSdfPath(smiles, idx, returnedPaths);
    return createViewerObject(smiles, sdfPath, nameMap);
  });
}

