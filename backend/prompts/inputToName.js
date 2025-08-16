// inputToName: build a prompt that converts arbitrary text input → canonical molecule names (optionally with PubChem CID)

function inputToName(objectDescription) {
  return `Task: Identify only fully specific chemical molecules present in the object. JSON response.

Object: ${objectDescription}

Requirements:
- Output ONLY fully specific molecule names (IUPAC or unique common names). No classes/families (e.g., "polyphenols", "sugars", "lipids"), no brand names.
- Prefer the exact isomer/form when typical for the object (e.g., "chlorogenic acid", "(+)-catechin", "L-ascorbic acid").
- If known, include PubChem CID (integer). If unknown, set cid to null.
- Include both major and minor constituents that are characteristic.
- Use canonical naming: prefer PubChem Title; if unavailable, use IUPAC; otherwise widely accepted common name.
- Exclude quantities, categories, or commentary.

Schema:
{
  "object": "string",
  "molecules": [
    {"name": "string", "cid": number|null}
  ]
}`;
}

module.exports = { inputToName };
