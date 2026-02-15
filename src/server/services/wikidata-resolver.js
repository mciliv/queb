const USER_AGENT = 'MolecularContents/1.0 (+github.com/example)';

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

async function fetchJson(url, headers = {}) {
  const f = await getFetch();
  const resp = await f(url, { headers: { 'User-Agent': USER_AGENT, ...headers } });
  if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
  return resp.json();
}

async function postSparql(query) {
  const f = await getFetch();
  const resp = await f('https://query.wikidata.org/sparql', {
    method: 'POST',
    headers: {
      'Content-Type': 'application/sparql-query',
      'Accept': 'application/sparql-results+json',
      'User-Agent': USER_AGENT
    },
    body: query
  });
  if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
  return resp.json();
}

async function searchEntity(title) {
  if (!title || typeof title !== 'string') return null;
  const enc = encodeURIComponent(title.trim());
  const url = `https://www.wikidata.org/w/api.php?action=wbsearchentities&search=${enc}&language=en&format=json`;
  const data = await fetchJson(url);
  const hit = (data?.search || [])[0];
  return hit?.id || null; // e.g., Q30021
}

async function fetchConstituents(qid) {
  if (!qid) return [];
  const query = `
SELECT ?part ?partLabel ?smiles ?inchikey WHERE {
  wd:${qid} wdt:P527 ?part .
  OPTIONAL { ?part wdt:P233 ?smiles. }
  OPTIONAL { ?part wdt:P235 ?inchikey. }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}`;
  const json = await postSparql(query);
  const rows = json?.results?.bindings || [];
  return rows.map(r => ({
    id: r.part?.value || null,
    label: r.partLabel?.value || null,
    smiles: r.smiles?.value || null,
    inchikey: r.inchikey?.value || null
  })).filter(x => x.label);
}

async function findObjectConstituents(objectTitle) {
  try {
    const qid = await searchEntity(objectTitle);
    if (!qid) return [];
    const parts = await fetchConstituents(qid);
    return parts;
  } catch (_) {
    return [];
  }
}

module.exports = {
  searchEntity,
  fetchConstituents,
  findObjectConstituents
};

