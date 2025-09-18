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
  const conversionMethods = ['openeye', 'rdkit', 'pubchem_cid', 'pubchem_name', 'llm'];
  
  for (const mol of input) {
    const name = typeof mol?.name === 'string' ? mol.name.trim() : '';
    let cid = mol?.cid ?? null;
    let smiles = null;
    let method = null;

    if (name) {
      // Try multiple methods with progressive fallback
      const methods = [
        { name: 'openeye', fn: () => convertNameToSmilesOpenEye(name) },
        { name: 'rdkit', fn: () => convertNameToSmilesRDKit(name) }
      ];
      
      for (const methodInfo of methods) {
        if (!smiles) {
          try {
            const result = await methodInfo.fn();
            if (result) {
              smiles = result;
              method = methodInfo.name;
              break;
            }
          } catch (error) {
            console.warn(`${methodInfo.name} conversion failed for ${name}:`, error.message);
          }
        }
      }

      // PubChem fallback with error handling
      if (!smiles) {
        try {
          let resolveNameFn = null;
          let getPropertiesByCIDFn = null;
          try {
            const mod = require('../services/pubchem');
            resolveNameFn = typeof mod?.resolveName === 'function' ? mod.resolveName : null;
            getPropertiesByCIDFn = typeof mod?.getPropertiesByCID === 'function' ? mod.getPropertiesByCID : null;
          } catch (pubchemError) {
            console.warn(`PubChem module unavailable for ${name}:`, pubchemError.message);
          }

          if (cid && getPropertiesByCIDFn) {
            try {
              const props = await getPropertiesByCIDFn(cid);
              const found = props?.smiles || null;
              if (found) {
                smiles = found;
                method = 'pubchem_cid';
              }
            } catch (cidError) {
              console.warn(`PubChem CID lookup failed for ${name} (CID: ${cid}):`, cidError.message);
            }
          }
          
          if (!smiles && resolveNameFn) {
            try {
              const res = await resolveNameFn(name);
              if (res) {
                cid = res?.cid || cid;
                if (res?.smiles) {
                  smiles = res.smiles;
                  method = 'pubchem_name';
                }
              }
            } catch (nameError) {
              console.warn(`PubChem name lookup failed for ${name}:`, nameError.message);
            }
          }
        } catch (pubchemGeneralError) {
          console.warn(`PubChem lookup general error for ${name}:`, pubchemGeneralError.message);
        }
      }

      // LLM fallback with enhanced error handling
      if (!smiles && llmClient) {
        try {
          const prompt = buildNameToSmilesPrompt({ object, molecules: [{ name, cid }] });
          const response = await llmClient.chat.completions.create({
            model: process.env.OPENAI_MODEL || process.env.OPENAI_DEFAULT_MODEL || 'gpt-4o',
            messages: [{ role: 'user', content: prompt }],
            max_tokens: 500,
            temperature: 0.1,
            response_format: { type: 'json_object' },
          });
          
          if (response && response.choices?.[0]?.message?.content) {
            try {
              const parsed = JSON.parse(response.choices[0].message.content);
              const llmResult = parsed.molecules?.[0];
              if (llmResult?.smiles) {
                // Basic SMILES validation
                const smilesPattern = /^[A-Za-z0-9@+\-\[\]()=#\/\\]+$/;
                if (smilesPattern.test(llmResult.smiles)) {
                  smiles = llmResult.smiles;
                  method = 'llm';
                } else {
                  console.warn(`Invalid SMILES format from LLM for ${name}: ${llmResult.smiles}`);
                }
              }
            } catch (parseError) {
              console.warn(`LLM response parsing failed for ${name}:`, parseError.message);
            }
          }
        } catch (llmError) {
          console.warn(`LLM conversion failed for ${name}:`, llmError.message);
        }
      }
    }

    results.push({
      name,
      cid,
      smiles: smiles || null,
      status: smiles ? 'ok' : 'lookup_required',
      method: method || 'failed',
      attempted_methods: conversionMethods.filter(m => 
        (m === 'openeye' || m === 'rdkit') ? true :
        (m === 'pubchem_cid' && cid) ? true :
        (m === 'pubchem_name') ? true :
        (m === 'llm' && llmClient) ? true : false
      ),
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
