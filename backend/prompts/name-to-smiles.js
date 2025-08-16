// Name to SMILES conversion using algorithmic libraries with LLM backup

const { spawn } = require('child_process');

async function convertNameToSmilesRDKit(name) {
  return new Promise((resolve) => {
    const pythonScript = `
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import sys

name = "${name.replace(/"/g, '\\"')}"
try:
    mol = Chem.MolFromName(name)
    if mol:
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        print(smiles)
    else:
        print("NONE")
except:
    print("NONE")
`;
    
    const python = spawn('python3', ['-c', pythonScript]);
    let output = '';
    
    python.stdout.on('data', (data) => {
      output += data.toString();
    });
    
    python.on('close', (code) => {
      const smiles = output.trim();
      resolve(smiles === 'NONE' || !smiles ? null : smiles);
    });
    
    python.on('error', () => {
      resolve(null);
    });
  });
}

async function convertNameToSmilesOpenEye(name) {
  return new Promise((resolve) => {
    const pythonScript = `
try:
    from openeye import oechem
    mol = oechem.OEGraphMol()
    if oechem.OEParseIUPACName(mol, "${name.replace(/"/g, '\\"')}"):
        smiles = oechem.OEMolToSmiles(mol)
        print(smiles)
    else:
        print("NONE")
except:
    print("NONE")
`;
    
    const python = spawn('python3', ['-c', pythonScript]);
    let output = '';
    
    python.stdout.on('data', (data) => {
      output += data.toString();
    });
    
    python.on('close', (code) => {
      const smiles = output.trim();
      resolve(smiles === 'NONE' || !smiles ? null : smiles);
    });
    
    python.on('error', () => {
      resolve(null);
    });
  });
}

async function convertNamesToSmiles(payload, llmClient = null) {
  const object = payload?.object || "";
  const input = Array.isArray(payload?.molecules) ? payload.molecules : [];
  const results = [];
  
  for (const mol of input) {
    const name = mol?.name || '';
    let cid = mol?.cid ?? null;
    let smiles = null;
    let method = null;
    
    // 1. Try OpenEye toolkit first (most comprehensive)
    try {
      smiles = await convertNameToSmilesOpenEye(name);
      if (smiles) method = 'openeye';
    } catch (_) {}
    
    // 2. Try RDKit if OpenEye fails
    if (!smiles) {
      try {
        smiles = await convertNameToSmilesRDKit(name);
        if (smiles) method = 'rdkit';
      } catch (_) {}
    }
    
    // 3. Try PubChem programmatic lookup
    if (!smiles) {
      try {
        const { resolveName, getPropertiesByCID } = require('../services/pubchem');
        
        if (cid) {
          const props = await getPropertiesByCID(cid);
          smiles = props?.smiles || null;
          if (smiles) method = 'pubchem_cid';
        }
        
        if (!smiles) {
          const res = await resolveName(name);
          cid = res?.cid || cid;
          smiles = res?.smiles || null;
          if (smiles) method = 'pubchem_name';
        }
      } catch (_) {}
    }
    
    // 4. LLM backup for edge cases not in databases
    if (!smiles && llmClient) {
      try {
        const prompt = buildNameToSmilesPrompt({ object, molecules: [{ name, cid }] });
        const response = await llmClient.chat.completions.create({
          model: "gpt-4o",
          messages: [{ role: "user", content: prompt }],
          max_tokens: 500,
          temperature: 0.1,
          response_format: { type: "json_object" }
        });
        
        const parsed = JSON.parse(response.choices[0].message.content);
        const llmResult = parsed.molecules?.[0];
        if (llmResult?.smiles) {
          smiles = llmResult.smiles;
          method = 'llm';
        }
      } catch (_) {}
    }
    
    results.push({
      name,
      cid,
      smiles: smiles || null,
      status: smiles ? 'ok' : 'lookup_required',
      method: method || 'failed'
    });
  }
  
  return { object, molecules: results };
}

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

module.exports = { buildNameToSmilesPrompt, convertNamesToSmiles };
