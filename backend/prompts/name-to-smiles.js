// Name to SMILES conversion prompt builder

function buildNameToSmilesPrompt(payload) {
  const input = typeof payload === 'string' ? payload : JSON.stringify(payload, null, 2);
  return `Task: Convert each molecule to an isomeric SMILES. JSON response.

Rules:
- For entries with a PubChem CID, use the PubChem IsomericSMILES for that CID.
- Without CID, provide the best-known isomeric SMILES for the named compound.
- If uncertain or ambiguous, set "smiles": null and "status": "lookup_required".
- Return JSON only, no prose.

Input:
${input}

Output schema:
{
  "object": "string",
  "molecules": [
    {"name": "string", "cid": number|null, "smiles": "string|null", "status": "ok|lookup_required"}
  ]
}`;
}

module.exports = { buildNameToSmilesPrompt };


