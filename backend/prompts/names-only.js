// Names-only extraction prompt builder

function buildNamesOnlyPrompt(objectDescription) {
  return `Task: Extract specific molecule names only. JSON response.

Object: ${objectDescription}

Rules:
- Only fully specific molecules (no categories or classes)
- Prefer canonical PubChem Title; otherwise IUPAC/common
- Include PubChem CID if confidently known, else null

Schema:
{
  "object": "string",
  "molecules": [
    {"name": "string", "cid": number|null}
  ]
}`;
}

function listChemicalsFromTextSpecification(text) {
  return buildNamesOnlyPrompt(text);
}

module.exports = { buildNamesOnlyPrompt, listChemicalsFromTextSpecification };


