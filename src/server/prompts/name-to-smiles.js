// Name to SMILES conversion using external Python helpers with LLM backup

const { spawn } = require('child_process');
const path = require('path');
const logger = require('../services/file-logger') || console;

// Configuration with environment-specific defaults
const CONFIG = {
  PY_HELPER: path.join(__dirname, '..', 'python', 'name_to_smiles.py'),
  PY_TIMEOUT_MS: parseInt(process.env.NAME_TO_SMILES_TIMEOUT_MS) || 15000,
  MAX_RETRIES: parseInt(process.env.NAME_TO_SMILES_MAX_RETRIES) || 2,
  RETRY_DELAY_MS: parseInt(process.env.NAME_TO_SMILES_RETRY_DELAY_MS) || 1000,
  CACHE_TTL_MS: parseInt(process.env.NAME_TO_SMILES_CACHE_TTL_MS) || 3600000, // 1 hour
  MAX_CACHE_SIZE: parseInt(process.env.NAME_TO_SMILES_MAX_CACHE_SIZE) || 1000
};

// Enhanced caching with TTL and size limits
const conversionCache = new Map();
const cacheMetadata = new Map(); // stores { timestamp, accessCount }

function cleanupCache() {
  const now = Date.now();
  const entries = Array.from(cacheMetadata.entries());
  
  // Remove expired entries
  entries.forEach(([key, meta]) => {
    if (now - meta.timestamp > CONFIG.CACHE_TTL_MS) {
      conversionCache.delete(key);
      cacheMetadata.delete(key);
    }
  });
  
  // If still over size limit, remove least recently accessed
  if (conversionCache.size > CONFIG.MAX_CACHE_SIZE) {
    const sortedEntries = entries
      .filter(([key]) => conversionCache.has(key))
      .sort((a, b) => a[1].accessCount - b[1].accessCount)
      .slice(0, conversionCache.size - CONFIG.MAX_CACHE_SIZE);
    
    sortedEntries.forEach(([key]) => {
      conversionCache.delete(key);
      cacheMetadata.delete(key);
    });
  }
}

function getCachedResult(key) {
  if (conversionCache.has(key)) {
    const meta = cacheMetadata.get(key);
    if (meta && Date.now() - meta.timestamp < CONFIG.CACHE_TTL_MS) {
      meta.accessCount++;
      return conversionCache.get(key);
    }
    // Expired
    conversionCache.delete(key);
    cacheMetadata.delete(key);
  }
  return null;
}

function setCachedResult(key, result) {
  cleanupCache();
  conversionCache.set(key, result);
  cacheMetadata.set(key, { timestamp: Date.now(), accessCount: 1 });
}

// Custom error classes for better error handling
class NameToSmilesError extends Error {
  constructor(message, type, originalError = null) {
    super(message);
    this.name = 'NameToSmilesError';
    this.type = type;
    this.originalError = originalError;
  }
}

class ValidationError extends NameToSmilesError {
  constructor(message, originalError = null) {
    super(message, 'VALIDATION_ERROR', originalError);
  }
}

class TimeoutError extends NameToSmilesError {
  constructor(message, originalError = null) {
    super(message, 'TIMEOUT_ERROR', originalError);
  }
}

class PythonExecutionError extends NameToSmilesError {
  constructor(message, originalError = null) {
    super(message, 'PYTHON_EXECUTION_ERROR', originalError);
  }
}

class NetworkError extends NameToSmilesError {
  constructor(message, originalError = null) {
    super(message, 'NETWORK_ERROR', originalError);
  }
}

async function runPythonNameToSmiles(toolkit, rawName, timeoutMs = CONFIG.PY_TIMEOUT_MS) {
  // Input validation
  if (!rawName || typeof rawName !== 'string') {
    throw new ValidationError('Invalid input: name must be a non-empty string');
  }
  
  const safeName = rawName.trim();
  if (safeName.length === 0) {
    throw new ValidationError('Invalid input: name cannot be empty after trimming');
  }
  
  if (safeName.length > 1000) {
    throw new ValidationError('Invalid input: name too long (max 1000 characters)');
  }
  
  // Check cache first
  const cacheKey = `${toolkit}:${safeName.toLowerCase()}`;
  const cached = getCachedResult(cacheKey);
  if (cached !== null) {
    logger.debug(`Cache hit for ${toolkit} conversion of "${safeName}"`);
    return cached;
  }
  
  logger.info(`Converting "${safeName}" to SMILES using ${toolkit}`);
  
  return new Promise((resolve, reject) => {
    const args = [CONFIG.PY_HELPER, '--toolkit', toolkit, '--name', safeName];
    const python = spawn('python3', args, { 
      stdio: ['ignore', 'pipe', 'pipe'],
      timeout: timeoutMs 
    });
    
    let output = '';
    let errorOutput = '';
    let settled = false;

    const finish = (result, error = null) => {
      if (settled) return;
      settled = true;
      clearTimeout(timer);
      
      if (error) {
        logger.error(`Python ${toolkit} conversion failed for "${safeName}":`, error.message);
        reject(error);
      } else {
        // Cache the result
        setCachedResult(cacheKey, result);
        logger.debug(`Successfully converted "${safeName}" using ${toolkit}: ${result ? 'SUCCESS' : 'NO_RESULT'}`);
        resolve(result);
      }
    };

    const timer = setTimeout(() => {
      try { 
        python.kill('SIGKILL'); 
      } catch (killError) {
        logger.warn(`Failed to kill Python process: ${killError.message}`);
      }
      finish(null, new TimeoutError(`Python ${toolkit} conversion timed out after ${timeoutMs}ms for "${safeName}"`));
    }, timeoutMs);

    python.stdout.on('data', (data) => {
      output += data.toString();
    });
    
    python.stderr.on('data', (data) => {
      errorOutput += data.toString();
    });

    python.on('close', (code, signal) => {
      if (signal === 'SIGKILL') {
        return; // Timeout already handled
      }
      
      const smiles = (output || '').trim();
      
      if (code !== 0) {
        const errorMsg = errorOutput.trim() || `Python process exited with code ${code}`;
        finish(null, new PythonExecutionError(`${toolkit} conversion failed: ${errorMsg}`));
        return;
      }
      
      const result = (smiles === 'NONE' || smiles.length === 0) ? null : smiles;
      
      // Basic SMILES validation
      if (result && !isValidSmilesFormat(result)) {
        logger.warn(`Invalid SMILES format returned by ${toolkit} for "${safeName}": ${result}`);
        finish(null, new ValidationError(`Invalid SMILES format returned: ${result}`));
        return;
      }
      
      finish(result);
    });

    python.on('error', (err) => {
      finish(null, new PythonExecutionError(`Failed to spawn Python process: ${err.message}`, err));
    });
  });
}

function isValidSmilesFormat(smiles) {
  if (!smiles || typeof smiles !== 'string') return false;
  
  // Basic SMILES pattern validation
  const smilesPattern = /^[A-Za-z0-9@+\-\[\]()=#\/\\.*:]+$/;
  if (!smilesPattern.test(smiles)) return false;
  
  // Check for balanced brackets
  const brackets = smiles.match(/[\[\]()]/g) || [];
  let squareCount = 0, parenCount = 0;
  
  for (const bracket of brackets) {
    if (bracket === '[') squareCount++;
    else if (bracket === ']') squareCount--;
    else if (bracket === '(') parenCount++;
    else if (bracket === ')') parenCount--;
    
    if (squareCount < 0 || parenCount < 0) return false;
  }
  
  return squareCount === 0 && parenCount === 0;
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
 
  logger.info(`Starting conversion of ${input.length} molecules to SMILES`);
 
  for (const mol of input) {
    const name = typeof mol?.name === 'string' ? mol.name.trim() : '';
    let smiles = null;
    let method = null;
    const errors = [];
 
    if (name) {
      const attempts = [
        { label: 'openeye', fn: () => convertNameToSmilesOpenEye(name) },
        { label: 'rdkit', fn: () => convertNameToSmilesRDKit(name) }
      ];
 
      for (const attempt of attempts) {
        if (smiles) break;
        try {
          const result = await attempt.fn();
          if (result) {
            smiles = result;
            method = attempt.label;
            logger.info(`Successfully converted "${name}" using ${method}`);
          }
        } catch (err) {
          errors.push(`${attempt.label}: ${err.message}`);
          logger.warn(`${attempt.label} conversion failed for "${name}":`, err.message);
        }
      }
 
      if (!smiles) {
        const resolved = await resolveSmilesViaPubChem(name, errors);
        if (resolved) {
          smiles = resolved;
          method = 'pubchem_name';
          logger.info(`Successfully resolved "${name}" via PubChem name lookup`);
        }
      }
 
      if (!smiles && llmClient) {
        try {
          const prompt = buildNameToSmilesPrompt({ object, molecules: [{ name }] });
          const response = await llmClient.chat.completions.create({
            model: process.env.OPENAI_MODEL || process.env.OPENAI_DEFAULT_MODEL || 'gpt-4o',
            messages: [{ role: 'user', content: prompt }],
            max_tokens: 500,
            response_format: { type: 'json_object' },
          });
 
          const content = response?.choices?.[0]?.message?.content;
          if (content) {
            const parsed = JSON.parse(content);
            const llmSmiles = parsed?.molecules?.[0]?.smiles;
            if (llmSmiles && isValidSmilesFormat(llmSmiles)) {
              smiles = llmSmiles;
              method = 'llm';
              logger.info(`Successfully converted "${name}" using LLM`);
            } else if (llmSmiles) {
              errors.push(`llm: Invalid SMILES format returned: ${llmSmiles}`);
              logger.warn(`Invalid SMILES format from LLM for ${name}: ${llmSmiles}`);
            }
          }
        } catch (err) {
          errors.push(`llm: ${err.message}`);
          logger.warn(`LLM conversion failed for ${name}:`, err.message);
        }
      }
    }
 
    results.push({
      name,
      smiles: smiles || null,
      status: smiles ? 'ok' : 'lookup_required',
      method: method || null,
      errors: errors.length ? errors : undefined,
    });
  }
 
  const successCount = results.filter(r => r.smiles).length;
  logger.info(`Conversion complete: ${successCount}/${results.length} molecules successfully converted to SMILES`);
 
  return {
    object,
    molecules: results,
  };
}

async function resolveSmilesViaPubChem(name, errors) {
  try {
    const mod = require('../services/pubchem');
    if (typeof mod?.resolveName === 'function') {
      const res = await mod.resolveName(name);
      if (res?.smiles && isValidSmilesFormat(res.smiles)) {
        return res.smiles;
      }
    }
  } catch (err) {
    if (err.code !== 'MODULE_NOT_FOUND') {
      errors.push(`pubchem_module: ${err.message}`);
      logger.warn(`PubChem module lookup failed for ${name}:`, err.message);
    }
  }
 
  try {
    const mod = require('../services/name-resolver');
    if (typeof mod?.resolveName === 'function') {
      const res = await mod.resolveName(name);
      if (res?.smiles && isValidSmilesFormat(res.smiles)) {
        return res.smiles;
      }
    }
  } catch (err) {
    errors.push(`pubchem_name: ${err.message}`);
    logger.warn(`Name resolver lookup failed for ${name}:`, err.message);
  }
 
  return null;
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
    {"name": "string", "smiles": "string|null", "status": "ok|lookup_required"}
  ]
}`;
}

module.exports = { buildNameToSmilesPrompt, convertNamesToSmiles };
