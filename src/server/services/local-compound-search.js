/**
 * local-compound-search.js - Local compound search service
 *
 * This service provides local search capabilities for chemical compounds based on object names.
 * It maintains a local database/cache of compounds that can be searched without external API calls,
 * providing fuzzy matching, ranking, and fast retrieval of compound information.
 *
 * Key features:
 * - Local compound database with fuzzy search
 * - Multiple data sources (PubChem, ChEMBL, local files)
 * - Search result ranking and filtering
 * - Persistent storage with JSON file cache
 * - Integration with existing name resolver
 */

const fs = require('fs');
const path = require('path');
const crypto = require('crypto');

class LocalCompoundSearch {
  constructor(options = {}) {
    this.dataDir = options.dataDir || path.join(__dirname, '..', '..', '..', 'data');
    this.compoundsFile = path.join(this.dataDir, 'local-compounds.json');
    this.indexFile = path.join(this.dataDir, 'search-index.json');

    // Search configuration
    this.searchConfig = {
      maxResults: options.maxResults || 50,
      minScore: options.minScore || 0.2,
      fuzzyThreshold: options.fuzzyThreshold || 0.6,
      enableCaching: options.enableCaching !== false,
      ...options
    };

    // In-memory data structures
    this.compounds = new Map(); // cid -> compound data
    this.searchIndex = new Map(); // search term -> [cid, score] pairs
    this.nameIndex = new Map(); // name variations -> cid

    // Ensure data directory exists
    this.ensureDataDirectory();

    // Load existing data
    if (this.searchConfig.enableCaching) {
      this.loadLocalData();
    }
  }

  ensureDataDirectory() {
    if (!fs.existsSync(this.dataDir)) {
      fs.mkdirSync(this.dataDir, { recursive: true });
    }
  }

  /**
   * Load compound data and search index from local files
   */
  loadLocalData() {
    try {
      // Load compounds database
      if (fs.existsSync(this.compoundsFile)) {
        const compoundsData = JSON.parse(fs.readFileSync(this.compoundsFile, 'utf8'));
        this.compounds = new Map(compoundsData);
        console.log(`Loaded ${this.compounds.size} compounds from local cache`);
      }

      // Load search index
      if (fs.existsSync(this.indexFile)) {
        const indexData = JSON.parse(fs.readFileSync(this.indexFile, 'utf8'));
        this.searchIndex = new Map(indexData);
        console.log(`Loaded search index with ${this.searchIndex.size} terms`);
      }
    } catch (error) {
      console.warn('Failed to load local compound data:', error.message);
    }
  }

  /**
   * Save compound data and search index to local files
   */
  saveLocalData() {
    if (!this.searchConfig.enableCaching) return;

    try {
      // Save compounds database
      const compoundsData = Array.from(this.compounds.entries());
      fs.writeFileSync(this.compoundsFile, JSON.stringify(compoundsData, null, 2));

      // Save search index
      const indexData = Array.from(this.searchIndex.entries());
      fs.writeFileSync(this.indexFile, JSON.stringify(indexData, null, 2));

      console.log(`Saved ${this.compounds.size} compounds and ${this.searchIndex.size} search terms`);
    } catch (error) {
      console.error('Failed to save local compound data:', error.message);
    }
  }

  /**
   * Add a compound to the local database
   * @param {Object} compound - Compound data
   * @param {string} compound.cid - PubChem CID
   * @param {string} compound.name - Primary name
   * @param {string[]} compound.names - Alternative names
   * @param {string} compound.smiles - SMILES notation
   * @param {string} compound.iupac - IUPAC name
   * @param {string} compound.formula - Molecular formula
   * @param {number} compound.molecularWeight - Molecular weight
   */
  addCompound(compound) {
    if (!compound.cid) {
      throw new Error('Compound must have a CID');
    }

    // Normalize compound data
    const normalizedCompound = {
      ...compound, // Spread first so original values can be overridden
      cid: compound.cid.toString(),
      name: compound.name || compound.title || '',
      names: Array.isArray(compound.names) ? compound.names : [],
      smiles: compound.smiles || '',
      iupac: compound.iupac || '',
      formula: compound.formula || '',
      molecularWeight: compound.molecularWeight || 0,
      source: compound.source || 'unknown',
      addedAt: compound.addedAt || new Date().toISOString()
    };

    // Add to compounds map
    this.compounds.set(normalizedCompound.cid, normalizedCompound);

    // Update search index
    this.updateSearchIndex(normalizedCompound);

    return normalizedCompound;
  }

  /**
   * Update the search index for a compound
   * @param {Object} compound - Normalized compound data
   */
  updateSearchIndex(compound) {
    const searchTerms = this.generateSearchTerms(compound);

    searchTerms.forEach(({ term, score }) => {
      const termKey = term.toLowerCase();
      if (!this.searchIndex.has(termKey)) {
        this.searchIndex.set(termKey, []);
      }

      const entries = this.searchIndex.get(termKey);
      const existingIndex = entries.findIndex(([cid]) => cid === compound.cid);

      if (existingIndex >= 0) {
        entries[existingIndex] = [compound.cid, score || 1.0];
      } else {
        entries.push([compound.cid, score || 1.0]);
      }

      // Sort by score (highest first)
      entries.sort((a, b) => b[1] - a[1]);
    });
  }

  /**
   * Generate search terms for a compound
   * @param {Object} compound - Compound data
   * @returns {Array} Array of search terms with scores
   */
  generateSearchTerms(compound) {
    const terms = [];

    // Primary name (highest score)
    if (compound.name) {
      terms.push({ term: compound.name, score: 1.0 });
    }

    // IUPAC name (high score)
    if (compound.iupac && compound.iupac !== compound.name) {
      terms.push({ term: compound.iupac, score: 0.9 });
    }

    // Alternative names (medium-high score)
    compound.names.forEach((name, index) => {
      if (name && name !== compound.name && name !== compound.iupac) {
        const score = Math.max(0.8 - (index * 0.05), 0.6); // Higher base score for alternative names
        terms.push({ term: name, score });

        // Also add individual words from alternative names
        const words = name.split(/\s+/);
        words.forEach(word => {
          if (word.length > 3) { // Require longer words for better matches
            terms.push({ term: word, score: score * 0.9 });
          }
        });
      }
    });

    // Chemical formula (medium score)
    if (compound.formula) {
      terms.push({ term: compound.formula, score: 0.6 });
    }

    // Generate partial terms for fuzzy matching
    terms.forEach(({ term, score }) => {
      // Split into words for partial matching
      const words = term.split(/\s+/);
      words.forEach(word => {
        if (word.length > 3) {
          terms.push({ term: word, score: score * 0.8 });
        }
      });

      // Add common abbreviations or variations
      if (term.includes('acid')) {
        terms.push({ term: term.replace('acid', '').trim(), score: score * 0.7 });
      }
      if (term.includes('chloride')) {
        terms.push({ term: term.replace('chloride', '').trim(), score: score * 0.7 });
      }
    });

    return terms;
  }

  /**
   * Search for compounds by name or query
   * @param {string} query - Search query
   * @param {Object} options - Search options
   * @returns {Array} Array of matching compounds with scores
   */
  search(query, options = {}) {
    if (!query || typeof query !== 'string') {
      return [];
    }

    const config = { ...this.searchConfig, ...options };
    const normalizedQuery = query.toLowerCase().trim();

    if (!normalizedQuery) {
      return [];
    }

    const matches = new Map(); // cid -> { compound, score, matchType }

    // Direct exact matches (highest priority)
    if (this.searchIndex.has(normalizedQuery)) {
      const entries = this.searchIndex.get(normalizedQuery);
      entries.forEach(([cid, score]) => {
        const compound = this.compounds.get(cid);
        if (compound) {
          matches.set(cid, {
            compound,
            score: score * 1.0, // Exact match bonus
            matchType: 'exact'
          });
        }
      });
    }

    // Fuzzy matching for partial terms
    this.searchIndex.forEach((entries, term) => {
      const similarity = this.calculateSimilarity(normalizedQuery, term);
      if (similarity >= config.fuzzyThreshold) {
        entries.forEach(([cid, termScore]) => {
          const compound = this.compounds.get(cid);
          if (compound) {
            const existing = matches.get(cid);
            const newScore = similarity * termScore * 0.9; // Fuzzy match penalty

            if (!existing || newScore > existing.score) {
              matches.set(cid, {
                compound,
                score: newScore,
                matchType: 'fuzzy'
              });
            }
          }
        });
      }
    });

    // Word-based matching for multi-word queries
    const queryWords = normalizedQuery.split(/\s+/);
    if (queryWords.length > 1) {
      queryWords.forEach(word => {
        if (this.searchIndex.has(word)) {
          const entries = this.searchIndex.get(word);
          entries.forEach(([cid, termScore]) => {
            const compound = this.compounds.get(cid);
            if (compound) {
              const existing = matches.get(cid);
              const newScore = termScore * 0.7; // Word match penalty

              if (!existing || newScore > existing.score) {
                matches.set(cid, {
                  compound,
                  score: newScore,
                  matchType: 'word'
                });
              }
            }
          });
        }
      });
    }

    // Filter by minimum score and convert to array
    let results = Array.from(matches.values())
      .filter(match => match.score >= config.minScore)
      .sort((a, b) => b.score - a.score)
      .slice(0, config.maxResults);

    return results;
  }

  /**
   * Calculate string similarity using Levenshtein distance
   * @param {string} str1 - First string
   * @param {string} str2 - Second string
   * @returns {number} Similarity score (0-1, 1 being identical)
   */
  calculateSimilarity(str1, str2) {
    if (str1 === str2) return 1.0;

    const len1 = str1.length;
    const len2 = str2.length;

    if (len1 === 0 || len2 === 0) return 0.0;

    // Simple length-based filtering for performance
    const maxLen = Math.max(len1, len2);
    const minLen = Math.min(len1, len2);
    if (maxLen > minLen * 2) return 0.0; // Too different in length

    // Levenshtein distance calculation
    const matrix = Array(len1 + 1).fill(null).map(() => Array(len2 + 1).fill(null));

    for (let i = 0; i <= len1; i++) matrix[i][0] = i;
    for (let j = 0; j <= len2; j++) matrix[0][j] = j;

    for (let i = 1; i <= len1; i++) {
      for (let j = 1; j <= len2; j++) {
        const cost = str1[i - 1] === str2[j - 1] ? 0 : 1;
        matrix[i][j] = Math.min(
          matrix[i - 1][j] + 1,     // deletion
          matrix[i][j - 1] + 1,     // insertion
          matrix[i - 1][j - 1] + cost // substitution
        );
      }
    }

    const distance = matrix[len1][len2];
    const maxDistance = Math.max(len1, len2);

    return 1 - (distance / maxDistance);
  }

  /**
   * Get compound by CID
   * @param {string|number} cid - PubChem CID
   * @returns {Object|null} Compound data or null if not found
   */
  getCompoundByCID(cid) {
    return this.compounds.get(cid.toString()) || null;
  }

  /**
   * Get all compounds (paginated)
   * @param {number} limit - Maximum number of results
   * @param {number} offset - Offset for pagination
   * @returns {Array} Array of compounds
   */
  getAllCompounds(limit = 100, offset = 0) {
    const allCompounds = Array.from(this.compounds.values());
    return allCompounds.slice(offset, offset + limit);
  }

  /**
   * Populate local database from PubChem API
   * @param {Array} compoundNames - Array of compound names to fetch
   * @returns {Promise<Object>} Results summary
   */
  async populateFromPubChem(compoundNames) {
    const results = {
      success: 0,
      failed: 0,
      skipped: 0,
      errors: []
    };

    // Import name resolver for PubChem access
    const { resolveName } = require('./name-resolver');

    for (const name of compoundNames) {
      try {
        // Check if we already have this compound
        const existing = Array.from(this.compounds.values())
          .find(c => c.name.toLowerCase() === name.toLowerCase() ||
                    c.names.some(n => n.toLowerCase() === name.toLowerCase()));

        if (existing) {
          results.skipped++;
          continue;
        }

        // Resolve compound from PubChem
        const resolved = await resolveName(name);

        if (resolved.smiles) {
          // Add to local database
          this.addCompound({
            cid: resolved.cid,
            name: resolved.title || name,
            names: [name], // Add search name as alternative
            smiles: resolved.smiles,
            iupac: resolved.iupac,
            source: 'pubchem'
          });

          results.success++;
        } else {
          results.failed++;
          results.errors.push(`No data found for: ${name}`);
        }

        // Small delay to be respectful to PubChem API
        await new Promise(resolve => setTimeout(resolve, 100));

      } catch (error) {
        results.failed++;
        results.errors.push(`${name}: ${error.message}`);
      }
    }

    // Save updated data
    this.saveLocalData();

    return results;
  }

  /**
   * Populate from a JSON file
   * @param {string} filePath - Path to JSON file containing compound data
   * @returns {Object} Results summary
   */
  populateFromFile(filePath) {
    const results = {
      success: 0,
      failed: 0,
      errors: []
    };

    try {
      const data = JSON.parse(fs.readFileSync(filePath, 'utf8'));
      const compounds = Array.isArray(data) ? data : [data];

      compounds.forEach(compound => {
        try {
          this.addCompound({
            ...compound,
            source: compound.source || 'file'
          });
          results.success++;
        } catch (error) {
          results.failed++;
          results.errors.push(`${compound.name || compound.cid}: ${error.message}`);
        }
      });

      this.saveLocalData();

    } catch (error) {
      results.failed++;
      results.errors.push(`File error: ${error.message}`);
    }

    return results;
  }

  /**
   * Get search statistics
   * @returns {Object} Database statistics
   */
  getStats() {
    return {
      totalCompounds: this.compounds.size,
      totalSearchTerms: this.searchIndex.size,
      dataSize: {
        compounds: this.compoundsFile,
        index: this.indexFile
      },
      cacheEnabled: this.searchConfig.enableCaching
    };
  }

  /**
   * Clear all local data
   */
  clear() {
    this.compounds.clear();
    this.searchIndex.clear();
    this.nameIndex.clear();

    // Remove files
    try {
      if (fs.existsSync(this.compoundsFile)) {
        fs.unlinkSync(this.compoundsFile);
      }
      if (fs.existsSync(this.indexFile)) {
        fs.unlinkSync(this.indexFile);
      }
    } catch (error) {
      console.warn('Failed to remove local data files:', error.message);
    }
  }

  /**
   * Rebuild search index from current compounds
   */
  rebuildIndex() {
    this.searchIndex.clear();
    this.nameIndex.clear();

    for (const compound of this.compounds.values()) {
      this.updateSearchIndex(compound);
    }

    this.saveLocalData();
    console.log(`Rebuilt search index for ${this.compounds.size} compounds`);
  }
}

module.exports = LocalCompoundSearch;
