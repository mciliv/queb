const fs = require('fs');
const path = require('path');
let toml; try { toml = require('toml'); } catch (_) { toml = null; }

function getConfig() {
  if (!toml) return {};
  const configPath = path.join(process.cwd(), 'config/config.toml');
  try {
    if (fs.existsSync(configPath)) {
      return toml.parse(fs.readFileSync(configPath, 'utf8'));
    }
  } catch (_) {}
  return {};
}

module.exports = { getConfig };






