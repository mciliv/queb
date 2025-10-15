const PUBCHEM_BASE = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug';
<<<<<<< Updated upstream
const OPSIN_BASE = 'https://opsin.ch.cam.ac.uk/opsin';
const CACTUS_BASE = 'https://cactus.nci.nih.gov/chemical/structure';
const CHEMBL_BASE = 'https://www.ebi.ac.uk/chembl/api/data';
=======
const logger = console;
>>>>>>> Stashed changes

// Enhanced caching with TTL and size limits
const CACHE_CONFIG = {
  TTL_MS: parseInt(process.env.PUBCHEM_CACHE_TTL_MS) || 1800000, // 30 minutes (reduced from 1 hour)
  MAX_SIZE: parseInt(process.env.PUBCHEM_MAX_CACHE_SIZE) || 500, // Reduced from 1000
  CLEANUP_INTERVAL_MS: parseInt(process.env.PUBCHEM_CACHE_CLEANUP_INTERVAL_MS) || 120000, // 2 minutes (reduced from 5 minutes)
  MAX_MEMORY_MB: parseInt(process.env.PUBCHEM_MAX_CACHE_MEMORY_MB) || 50 // 50MB memory limit for cache
};

// Cache with metadata for TTL and LRU
const namePropsCache = new Map(); // key: lowercased name -> { value: { smiles, title, iupac }, timestamp: number, accessCount: number }

// Periodic cache cleanup
let cleanupInterval = null;

// Start cleanup interval (can be stopped for testing)
function startCleanupInterval() {
  if (!cleanupInterval) {
    cleanupInterval = setInterval(() => {
      cleanupCache(namePropsCache);
    }, CACHE_CONFIG.CLEANUP_INTERVAL_MS);
    
    // Allow unref in Node.js to prevent keeping process alive
    if (cleanupInterval.unref) {
      cleanupInterval.unref();
    }
  }
}

// Stop cleanup interval (for testing)
function stopCleanupInterval() {
  if (cleanupInterval) {
    clearInterval(cleanupInterval);
    cleanupInterval = null;
  }
}

// Auto-start in non-test environments
if (process.env.NODE_ENV !== 'test') {
  startCleanupInterval();
}

function cleanupCache(cache) {
  const now = Date.now();
  const entries = Array.from(cache.entries());
  
  // Remove expired entries
  entries.forEach(([key, data]) => {
    if (now - data.timestamp > CACHE_CONFIG.TTL_MS) {
      cache.delete(key);
    }
  });
  
  // Estimate cache memory usage (rough approximation)
  const estimatedMemoryMB = cache.size * 0.1; // Rough estimate: 100KB per entry
  
  // If over memory limit, remove least recently accessed
  if (estimatedMemoryMB > CACHE_CONFIG.MAX_MEMORY_MB) {
    const sortedEntries = entries
      .filter(([key]) => cache.has(key))
      .sort((a, b) => a[1].accessCount - b[1].accessCount);
    
    // Remove 25% of entries to free up memory
    const toRemove = Math.floor(sortedEntries.length * 0.25);
    sortedEntries.slice(0, toRemove).forEach(([key]) => {
      cache.delete(key);
    });
  }
  
  // If still over size limit, remove least recently accessed
  if (cache.size > CACHE_CONFIG.MAX_SIZE) {
    const sortedEntries = entries
      .filter(([key]) => cache.has(key))
      .sort((a, b) => a[1].accessCount - b[1].accessCount)
      .slice(0, cache.size - CACHE_CONFIG.MAX_SIZE);
    
    sortedEntries.forEach(([key]) => {
      cache.delete(key);
    });
  }
}

function getCacheValue(cache, key) {
  const data = cache.get(key);
  if (data && Date.now() - data.timestamp < CACHE_CONFIG.TTL_MS) {
    data.accessCount++;
    return data.value;
  }
  if (data) {
    cache.delete(key); // Expired
  }
  return null;
}

function setCacheValue(cache, key, value) {
  cache.set(key, { value, timestamp: Date.now(), accessCount: 1 });
}

// Custom error classes
class PubChemError extends Error {
  constructor(message, statusCode = null, originalError = null) {
    super(message);
    this.name = 'PubChemError';
    this.statusCode = statusCode;
    this.originalError = originalError;
  }
}

class PubChemTimeoutError extends PubChemError {
  constructor(message, originalError = null) {
    super(message, 408, originalError);
    this.name = 'PubChemTimeoutError';
  }
}

class PubChemRateLimitError extends PubChemError {
  constructor(message, originalError = null) {
    super(message, 429, originalError);
    this.name = 'PubChemRateLimitError';
  }
}

// Enhanced retry/backoff configuration
const RETRY_CONFIG = {
  MAX_RETRIES: parseInt(process.env.PUBCHEM_MAX_RETRIES) || 3,
  BASE_DELAY_MS: parseInt(process.env.PUBCHEM_BASE_DELAY_MS) || 500,
  MAX_DELAY_MS: parseInt(process.env.PUBCHEM_MAX_DELAY_MS) || 5000,
  BACKOFF_MULTIPLIER: parseFloat(process.env.PUBCHEM_BACKOFF_MULTIPLIER) || 2.0
};

const sleep = (ms) => new Promise(r => setTimeout(r, ms));

function calculateDelay(attempt) {
  const delay = RETRY_CONFIG.BASE_DELAY_MS * Math.pow(RETRY_CONFIG.BACKOFF_MULTIPLIER, attempt - 1);
  return Math.min(delay, RETRY_CONFIG.MAX_DELAY_MS);
}

let cachedFetch = null;
async function getFetch() {
  if (cachedFetch) return cachedFetch;
  if (typeof fetch === 'function') {
    cachedFetch = (...args) => global.fetch(...args);
    return cachedFetch;
  }
  const mod = await import('node-fetch');
  cachedFetch = mod.default;
  return cachedFetch;
}

async function fetchJson(url, attempt = 1) {
  const f = await getFetch();
  
  for (let currentAttempt = attempt; currentAttempt <= RETRY_CONFIG.MAX_RETRIES; currentAttempt++) {
    try {
      logger.debug(`Fetching PubChem data from ${url} (attempt ${currentAttempt}/${RETRY_CONFIG.MAX_RETRIES})`);
      
      const resp = await f(url, {
        timeout: 10000,
        headers: {
          'User-Agent': 'MolecularAnalysisService/1.0',
          'Accept': 'application/json'
        }
      });
      
      if (!resp.ok) {
        const status = resp.status;
        const statusText = resp.statusText || 'Unknown Error';
        
        // Handle specific HTTP errors
        if (status === 404) {
          throw new PubChemError(`Resource not found: ${statusText}`, status);
        } else if (status === 429) {
          if (currentAttempt < RETRY_CONFIG.MAX_RETRIES) {
            const delay = calculateDelay(currentAttempt) * 2; // Extra delay for rate limiting
            logger.warn(`PubChem rate limit hit, retrying after ${delay}ms`);
            await sleep(delay);
            continue;
          }
          throw new PubChemRateLimitError(`Rate limit exceeded: ${statusText}`, status);
        } else if ([500, 502, 503, 504].includes(status)) {
          if (currentAttempt < RETRY_CONFIG.MAX_RETRIES) {
            const delay = calculateDelay(currentAttempt);
            logger.warn(`PubChem server error ${status}, retrying after ${delay}ms`);
            await sleep(delay);
            continue;
          }
          throw new PubChemError(`Server error: ${statusText}`, status);
        } else {
          throw new PubChemError(`HTTP error: ${statusText}`, status);
        }
      }
      
      const data = await resp.json();
      logger.debug(`Successfully fetched PubChem data from ${url}`);
      return data;
      
    } catch (err) {
      // Handle network and timeout errors
      if (err.name === 'AbortError' || err.code === 'ETIMEDOUT') {
        if (currentAttempt < RETRY_CONFIG.MAX_RETRIES) {
          const delay = calculateDelay(currentAttempt);
          logger.warn(`PubChem request timeout, retrying after ${delay}ms`);
          await sleep(delay);
          continue;
        }
        throw new PubChemTimeoutError(`Request timeout after ${RETRY_CONFIG.MAX_RETRIES} attempts`, err);
      } else if (err.name === 'FetchError' || err.code === 'ECONNRESET' || err.code === 'ENOTFOUND') {
        if (currentAttempt < RETRY_CONFIG.MAX_RETRIES) {
          const delay = calculateDelay(currentAttempt);
          logger.warn(`PubChem network error, retrying after ${delay}ms: ${err.message}`);
          await sleep(delay);
          continue;
        }
        throw new PubChemError(`Network error after ${RETRY_CONFIG.MAX_RETRIES} attempts: ${err.message}`, null, err);
      } else if (err instanceof PubChemError) {
        throw err; // Re-throw our custom errors
      } else {
        if (currentAttempt < RETRY_CONFIG.MAX_RETRIES) {
          const delay = calculateDelay(currentAttempt);
          logger.warn(`PubChem unexpected error, retrying after ${delay}ms: ${err.message}`);
          await sleep(delay);
          continue;
        }
        throw new PubChemError(`Unexpected error after ${RETRY_CONFIG.MAX_RETRIES} attempts: ${err.message}`, null, err);
      }
    }
  }
}
// Plain text fetch with basic retry
async function fetchText(url, attempt = 0) {
  const f = await getFetch();
  try {
    const resp = await f(url);
    if (!resp.ok) {
      const status = resp.status;
      if ((status === 429 || status === 500 || status === 502 || status === 503 || status === 504) && attempt < MAX_RETRIES) {
        await sleep(BASE_DELAY_MS * Math.pow(2, attempt));
        return fetchText(url, attempt + 1);
      }
      throw new Error(`HTTP ${status}`);
    }
    return resp.text();
  } catch (err) {
    if (attempt < MAX_RETRIES && (err.name === 'FetchError' || err.code === 'ECONNRESET' || err.code === 'ETIMEDOUT')) {
      await sleep(BASE_DELAY_MS * Math.pow(2, attempt));
      return fetchText(url, attempt + 1);
    }
    throw err;
  }
}

async function resolveCIDBySmiles(smiles) {
  const encoded = encodeURIComponent(smiles);
  const url = `${PUBCHEM_BASE}/compound/smiles/${encoded}/cids/JSON`;
  const data = await fetchJson(url).catch(() => null);
  const cids = data?.IdentifierList?.CID || [];
  return cids.length > 0 ? cids[0] : null;
}

async function resolveViaOPSIN(name) {
  try {
    const encoded = encodeURIComponent(name);
    const url = `${OPSIN_BASE}/${encoded}.json`;
    const data = await fetchJson(url);
    const smiles = data?.smiles || null;
    return smiles ? { smiles } : null;
  } catch (_) {
    return null;
  }
}

async function resolveViaCACTUS(name) {
  try {
    const encoded = encodeURIComponent(name);
    const url = `${CACTUS_BASE}/${encoded}/smiles`;
    const text = await fetchText(url);
    const smiles = (text || '').trim();
    if (!smiles || /not found/i.test(smiles)) return null;
    return { smiles };
  } catch (_) {
    return null;
  }
}

async function resolveViaChEMBL(name) {
  try {
    const q = encodeURIComponent(name);
    const url = `${CHEMBL_BASE}/molecule/search.json?q=${q}&limit=1`;
    const data = await fetchJson(url);
    const mol = data?.molecules?.[0];
    const smiles = mol?.molecule_structures?.canonical_smiles || null;
    return smiles ? { smiles, chembl_id: mol?.molecule_chembl_id || null, pref_name: mol?.pref_name || null } : null;
  } catch (_) {
    return null;
  }
}


async function resolveName(name) {
  // Returns { smiles, title, iupac } best-effort
  if (!name || typeof name !== 'string') {
    throw new PubChemError('Invalid name: must be a non-empty string');
  }
  
  const trimmedName = name.trim();
  if (trimmedName.length === 0) {
    throw new PubChemError('Invalid name: cannot be empty after trimming');
  }
  
  const cacheKey = trimmedName.toLowerCase();
  
  // Check cache first
  const cached = getCacheValue(namePropsCache, cacheKey);
  if (cached !== null) {
    logger.debug(`Cache hit for name resolution: "${trimmedName}"`);
    return cached;
  }
  
  logger.info(`Resolving name properties: "${trimmedName}"`);
  
  let smiles = null;
  let title = null;
  let iupac = null;
  const errors = [];
  
  // First try: direct properties lookup by name
  try {
    logger.debug(`Trying direct properties lookup for "${trimmedName}"`);
    const encoded = encodeURIComponent(trimmedName);
    const url = `${PUBCHEM_BASE}/compound/name/${encoded}/property/IsomericSMILES,CID,IUPACName,Title/JSON`;
    const data = await fetchJson(url);
    const props = data?.PropertyTable?.Properties?.[0] || {};
    smiles = props.IsomericSMILES || null;
    title = props.Title || null;
    iupac = props.IUPACName || null;
    
    if (smiles || title) {
      logger.info(`Direct properties lookup succeeded for "${trimmedName}"`);
    }
  } catch (directError) {
    errors.push(`Direct properties lookup failed: ${directError.message}`);
    logger.warn(`Direct properties lookup failed for "${trimmedName}":`, directError.message);
  }

  const result = { smiles, title, iupac, errors: errors.length > 0 ? errors : undefined };
  setCacheValue(namePropsCache, cacheKey, result);
  
  if (smiles) {
    logger.info(`Successfully resolved "${trimmedName}": SMILES=${smiles}`);
  } else {
    logger.warn(`Failed to resolve "${trimmedName}" to SMILES. Errors: ${errors.join('; ')}`);
  }
  
  return result;
}

async function downloadSDFBySmiles(smiles) {
  if (!smiles || typeof smiles !== 'string') {
    throw new Error('Invalid SMILES input');
  }
  
  const encoded = encodeURIComponent(smiles.trim());
  const url = `${PUBCHEM_BASE}/compound/SMILES/${encoded}/SDF`;
  
  try {
    const f = await getFetch();
    const resp = await f(url, { timeout: 10000 });
    
    if (!resp.ok) {
      if (resp.status === 404) {
        throw new Error(`Compound not found in PubChem`);
      }
      throw new Error(`PubChem API error: ${resp.status}`);
    }
    
    const text = await resp.text();
    if (!text || text.length < 100) {
      throw new Error('Invalid SDF response from PubChem');
    }
    
    return text;
  } catch (error) {
    if (error.name === 'AbortError' || error.code === 'ETIMEDOUT') {
      throw new Error('PubChem request timeout');
    }
    throw error;
  }
}

<<<<<<< Updated upstream
async function resolveName(name) {
  // Returns { cid, smiles, title, iupac } best-effort
  if (!name || typeof name !== 'string') {
    return { cid: null, smiles: null, title: null, iupac: null };
  }
  
  let cid = null;
  let title = null;
  let iupac = null;
  let smiles = null;
  const cacheKey = name.toLowerCase();
  if (namePropsCache.has(cacheKey)) return namePropsCache.get(cacheKey);
  
  try {
    cid = await resolveNameToCID(name);
  } catch (_) {}

  if (cid) {
    try {
      const p = await getPropertiesByCID(cid);
      smiles = p.smiles || null;
      title = p.title || null;
      iupac = p.iupac || null;
    } catch (_) {}
  } else {
    // Try direct properties by name
    try {
      const encoded = encodeURIComponent(name);
      const url = `${PUBCHEM_BASE}/compound/name/${encoded}/property/IsomericSMILES,CID,IUPACName,Title/JSON`;
      const data = await fetchJson(url);
      const props = data?.PropertyTable?.Properties?.[0] || {};
      cid = props.CID || null;
      smiles = props.IsomericSMILES || null;
      title = props.Title || null;
      iupac = props.IUPACName || null;
    } catch (_) {}
  }

  // If still unresolved, try alternative resolvers depending on configuration
  try {
    const configuration = require('../../core/Configuration');
    const chemCfg = configuration.get('chem') || { primary: 'pubchem', enableAlternates: false };
    if (!smiles && chemCfg.enableAlternates) {
      const altOpsin = await resolveViaOPSIN(name);
      if (altOpsin && altOpsin.smiles) {
        smiles = altOpsin.smiles;
      }
    }
    if (!smiles && chemCfg.enableAlternates) {
      const altCactus = await resolveViaCACTUS(name);
      if (altCactus && altCactus.smiles) {
        smiles = altCactus.smiles;
      }
    }
    if (!smiles && (chemCfg.primary === 'chembl' || chemCfg.enableAlternates)) {
      const altChEMBL = await resolveViaChEMBL(name);
      if (altChEMBL && altChEMBL.smiles) {
        smiles = altChEMBL.smiles;
        title = altChEMBL.pref_name || title;
      }
    }
  } catch (_) {}

  // If we obtained SMILES from alternates but no CID, try back-resolving CID via SMILES in PubChem
  if (smiles && !cid) {
    try {
      cid = await resolveCIDBySmiles(smiles);
      if (cid) {
        const p = await getPropertiesByCID(cid).catch(() => null);
        title = (p?.title) || title;
        iupac = (p?.iupac) || iupac;
      }
    } catch (_) {}
  }

  const result = { cid, smiles, title, iupac };
  namePropsCache.set(cacheKey, result);
  return result;
}

=======
>>>>>>> Stashed changes
module.exports = {
  resolveName,
  downloadSDFBySmiles,
<<<<<<< Updated upstream
  resolveCIDBySmiles,
  resolveViaOPSIN,
  resolveViaCACTUS,
  resolveViaChEMBL,
=======
  startCleanupInterval,
  stopCleanupInterval,
>>>>>>> Stashed changes
};


