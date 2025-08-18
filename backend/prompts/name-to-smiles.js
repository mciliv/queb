// Name to SMILES conversion using external Python helpers with LLM backup

const { spawn } = require('child_process');
const path = require('path');

const PY_HELPER = path.join(__dirname, '..', 'python', 'name_to_smiles.py');
const PY_TIMEOUT_MS = 15000;

function runPythonNameToSmiles(toolkit, rawName, timeoutMs = PY_TIMEOUT_MS) {
  const safeName = typeof rawName === 'string' ? rawName.trim() : '';
  if (safeName.length === 0) {
    return Promise.resolve(null);
  }
  return new Promise((resolve) => {
    const args = [PY_HELPER, '--toolkit', toolkit, '--name', safeName];
    const python = spawn('python3', args, { stdio: ['ignore', 'pipe', 'ignore'] });
    let output = '';
    let settled = false;

    const finish = (value) => {
      if (settled) return;
      settled = true;
      clearTimeout(timer);
      resolve(value);
    };

    const timer = setTimeout(() => {
      try { python.kill('SIGKILL'); } catch (_) {}
      finish(null);
    }, timeoutMs);

    python.stdout.on('data', (data) => {
      output += data.toString();
    });

    python.on('close', () => {
      const smiles = (output || '').trim();
      finish(smiles === 'NONE' || smiles.length === 0 ? null : smiles);
    });

    python.on('error', () => {
      finish(null);
    });
  });
}

async function convertNameToSmilesRDKit(name) {
  return runPythonNameToSmiles('rdkit', name);
}

async function convertNameToSmilesOpenEye(name) {
  return runPythonNameToSmiles('openeye', name);
}

async function convertNamesToSmiles(payload, llmClient = null) {
  const object = payload?.object || "";
  const input = Array.isArray(payload?.molecules) ? payload.molecules : [];
  const results = [];
  
  for (const mol of input) {
    const name = typeof mol?.name === 'string' ? mol.name.trim() : '';
    let cid = mol?.cid ?? null;
    let smiles = null;
    let method = null;

    if (name) {
      const openeye = await convertNameToSmilesOpenEye(name).catch(() => null);
      if (openeye) {
        smiles = openeye;
        method = 'openeye';
      }

      if (!smiles) {
        const rdkit = await convertNameToSmilesRDKit(name).catch(() => null);
        if (rdkit) {
          smiles = rdkit;
          method = 'rdkit';
        }
      }

      if (!smiles) {
        let resolveNameFn = null;
        let getPropertiesByCIDFn = null;
        try {
          const mod = require('../services/pubchem');
          resolveNameFn = typeof mod?.resolveName === 'function' ? mod.resolveName : null;
          getPropertiesByCIDFn = typeof mod?.getPropertiesByCID === 'function' ? mod.getPropertiesByCID : null;
        } catch (_) {}

        if (cid && getPropertiesByCIDFn) {
          const props = await getPropertiesByCIDFn(cid).catch(() => null);
          const found = props?.smiles || null;
          if (found) {
            smiles = found;
            method = 'pubchem_cid';
          }
        }
         if (!smiles && resolveNameFn) {
          const res = await resolveNameFn(name).catch(() => null);
          if (res) {
            cid = res?.cid || cid;
            if (res?.smiles) {
              smiles = res.smiles;
              method = 'pubchem_name';
            }
          }
        }
      }

      if (!smiles && llmClient) {
        const prompt = buildNameToSmilesPrompt({ object, molecules: [{ name, cid }] });
        const response = await llmClient.chat.completions
          .create({
            model: process.env.OPENAI_MODEL || process.env.OPENAI_DEFAULT_MODEL || 'gpt-4o',
            messages: [{ role: 'user', content: prompt }],
            max_tokens: 500,
            temperature: 0.1,
            response_format: { type: 'json_object' },
          })
          .catch(() => null);
        if (response && response.choices?.[0]?.message?.content) {
          try {
            const parsed = JSON.parse(response.choices[0].message.content);
            const llmResult = parsed.molecules?.[0];
            if (llmResult?.smiles) {
              smiles = llmResult.smiles;
              method = 'llm';
            }
          } catch (_) {}
        }
      }
    }

    results.push({
      name,
      cid,
      smiles: smiles || null,
      status: smiles ? 'ok' : 'lookup_required',
      method: method || 'failed',
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
