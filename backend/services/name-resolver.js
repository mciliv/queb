const fetch = require('node-fetch');

const PUBCHEM_BASE = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug';

async function fetchJson(url) {
  const resp = await fetch(url);
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
  const resp = await fetch(url);
  if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
  const text = await resp.text();
  return text; // caller writes to file
}

async function resolveName(name) {
  // Returns { cid, smiles, title, iupac } best-effort
  let cid = null;
  let title = null;
  let iupac = null;
  try {
    cid = await resolveNameToCID(name);
  } catch (_) {}

  let smiles = null;
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

  return { cid, smiles, title, iupac };
}

module.exports = {
  resolveName,
  resolveNameToCID,
  getPropertiesByCID,
  downloadSDFByCID,
};


