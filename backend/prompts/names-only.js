// Names-only extraction prompt builder

function buildNamesOnlyPrompt(objectDescription) {
  return `Task: Identify only fully specific chemical molecules present in the object. JSON response.

Object: ${objectDescription}

Requirements:
- Output ONLY fully specific molecule names (IUPAC or unique common names). No classes/families (e.g., "polyphenols", "sugars", "lipids"), no brand names.
- Prefer the exact isomer/form when typical for the object (e.g., "chlorogenic acid", "(+)-catechin", "L-ascorbic acid").
- If known, include PubChem CID (integer). If unknown, set cid to null.
- Include both major and minor constituents that are characteristic.
- Strictly no SMILES in this step.
- JSON response only, no prose.

Schema:
{
  "object": "string",
  "molecules": [
    {"name": "string", "cid": number|null}
  ]
}`;
}

module.exports = { buildNamesOnlyPrompt };


