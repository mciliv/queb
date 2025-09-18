function buildStructuralizePrompt(objectText) {
  return [
    'Task: Given an object description, list characteristic molecules with valid, verified SMILES. JSON only. Include a brief reason explaining your identification.',
    `Object: ${objectText}`,
    'Rules:\n- Use canonical names when possible\n- SMILES must be chemically plausible\n- If unknown, return an empty list',
    'JSON format:',
    '{ "object": "string", "reason": "string", "chemicals": [ { "name": "string", "smiles": "SMILES" } ] }'
  ].join('\n');
}

module.exports = { buildStructuralizePrompt };


