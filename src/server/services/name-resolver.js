const PUBCHEM_BASE = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug';

// Minimal in-memory caches to avoid hammering upstream on repeats
const nameToCIDCache = new Map(); // key: lowercased name -> cid|null
const cidPropsCache = new Map(); // key: cid -> { smiles, iupac, title }
const namePropsCache = new Map(); // key: lowercased name -> { cid, smiles, title, iupac }

// Simple retry/backoff for transient errors (e.g., 503)
const MAX_RETRIES = 2;
const BASE_DELAY_MS = 400;
const sleep = (ms) => new Promise(r => setTimeout(r, ms));

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

async function fetchJson(url, attempt = 0) {
  const f = await getFetch();
  try {
    const resp = await f(url);
    if (!resp.ok) {
      const status = resp.status;
      // Retry for transient upstream failures
      if ((status === 429 || status === 500 || status === 502 || status === 503 || status === 504) && attempt < MAX_RETRIES) {
        await sleep(BASE_DELAY_MS * Math.pow(2, attempt));
        return fetchJson(url, attempt + 1);
      }
      throw new Error(`HTTP ${status}`);
    }
    return resp.json();
  } catch (err) {
    // Retry on network errors
    if (attempt < MAX_RETRIES && (err.name === 'FetchError' || err.code === 'ECONNRESET' || err.code === 'ETIMEDOUT')) {
      await sleep(BASE_DELAY_MS * Math.pow(2, attempt));
      return fetchJson(url, attempt + 1);
    }
    throw err;
  }
}

async function resolveNameToCID(name) {
  const encoded = encodeURIComponent(name);
  const url = `${PUBCHEM_BASE}/compound/name/${encoded}/cids/JSON`;
  const key = String(name || '').toLowerCase();
  if (nameToCIDCache.has(key)) return nameToCIDCache.get(key);
  const data = await fetchJson(url);
  const cids = data?.IdentifierList?.CID || [];
  const cid = cids.length > 0 ? cids[0] : null;
  nameToCIDCache.set(key, cid);
  return cid;
}

async function getPropertiesByCID(cid) {
  const url = `${PUBCHEM_BASE}/compound/cid/${cid}/property/IsomericSMILES,IUPACName,Title/JSON`;
  if (cidPropsCache.has(cid)) return cidPropsCache.get(cid);
  const data = await fetchJson(url);
  const props = data?.PropertyTable?.Properties?.[0] || {};
  const out = {
    smiles: props.IsomericSMILES || null,
    iupac: props.IUPACName || null,
    title: props.Title || null
  };
  cidPropsCache.set(cid, out);
  return out;
}

async function downloadSDFByCID(cid) {
  const url = `${PUBCHEM_BASE}/compound/CID/${cid}/SDF`;
  const f = await getFetch();
  let attempt = 0;
  while (true) {
    const resp = await f(url);
    if (resp.ok) {
      const text = await resp.text();
      return text;
    }
    const status = resp.status;
    if (attempt < MAX_RETRIES && (status === 429 || status === 500 || status === 502 || status === 503 || status === 504)) {
      await sleep(BASE_DELAY_MS * Math.pow(2, attempt++));
      continue;
    }
    throw new Error(`HTTP ${status}`);
  }
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

  const result = { cid, smiles, title, iupac };
  namePropsCache.set(cacheKey, result);
  return result;
}

module.exports = {
  resolveName,
  resolveNameToCID,
  getPropertiesByCID,
  downloadSDFByCID,
  downloadSDFBySmiles,
};


