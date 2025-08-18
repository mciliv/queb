const PUBCHEM_BASE = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug';

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

async function fetchJson(url) {
  const f = await getFetch();
  const resp = await f(url);
  if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
  return resp.json();
}

async function resolveNameToCID(name) {
  const encoded = encodeURIComponent(name);
  const url = `${PUBCHEM_BASE}/compound/name/${encoded}/cids/JSON`;
  const data = await fetchJson(url);
  const cids = data?.IdentifierList?.CID || [];
  return cids.length > 0 ? cids[0] : null;
}

async function getPropertiesByCID(cid) {
  const url = `${PUBCHEM_BASE}/compound/cid/${cid}/property/IsomericSMILES,IUPACName,Title/JSON`;
  const data = await fetchJson(url);
  const props = data?.PropertyTable?.Properties?.[0] || {};
  return { 
    smiles: props.IsomericSMILES || null,
    iupac: props.IUPACName || null,
    title: props.Title || null
  };
}

async function downloadSDFByCID(cid) {
  const url = `${PUBCHEM_BASE}/compound/CID/${cid}/SDF`;
  const f = await getFetch();
  const resp = await f(url);
  if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
  const text = await resp.text();
  return text; // caller writes to file
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
  
  try {
    cid = await resolveNameToCID(name);
  } catch (error) {
    console.log(`CID lookup failed for ${name}: ${error.message}`);
  }

  if (cid) {
    try {
      const p = await getPropertiesByCID(cid);
      smiles = p.smiles || null;
      title = p.title || null;
      iupac = p.iupac || null;
    } catch (error) {
      console.log(`Properties lookup failed for CID ${cid}: ${error.message}`);
    }
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
    } catch (error) {
      console.log(`Direct name lookup failed for ${name}: ${error.message}`);
    }
  }

  return { cid, smiles, title, iupac };
}

module.exports = {
  resolveName,
  resolveNameToCID,
  getPropertiesByCID,
  downloadSDFByCID,
  downloadSDFBySmiles,
};


