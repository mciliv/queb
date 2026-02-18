const { execFile } = require('child_process');
const path = require('path');

const DEFAULT_TIMEOUT_MS = 10000;
const DEFAULT_LIMIT = 5;

function buildArgs({ city, lat, lng, limit, ingredientsPath, citiesPath }) {
  const args = [];

  if (city) {
    args.push('--city', city);
  } else if (lat !== undefined && lng !== undefined) {
    args.push('--lat', String(lat), '--lng', String(lng));
  } else {
    throw new Error('Provide city or lat/lng to fetch ingredients.');
  }

  if (limit) {
    args.push('--limit', String(limit));
  }

  if (ingredientsPath) {
    args.push('--ingredients', ingredientsPath);
  }

  if (citiesPath) {
    args.push('--cities', citiesPath);
  }

  return args;
}

function resolveRubyPaths() {
  const root = path.join(__dirname, '..', '..', '..');
  return {
    binary: path.join(root, 'ruby_app', 'bin', 'ingredient_delivery'),
    ingredientsPath: path.join(root, 'ruby_app', 'data', 'ingredients.json'),
    citiesPath: path.join(root, 'ruby_app', 'data', 'cities.json')
  };
}

function fetchHighestQualityIngredients({ city, lat, lng, limit = DEFAULT_LIMIT, timeoutMs = DEFAULT_TIMEOUT_MS }) {
  const { binary, ingredientsPath, citiesPath } = resolveRubyPaths();

  return new Promise((resolve, reject) => {
    let args;
    try {
      // Block: fetchHighestQualityIngredients.
      // Refined prompt: integrate a Ruby app into the existing pipeline so it can
      // deliver highest-quality ingredients for a user's location with minimal
      // latency and a clear request contract. Prompt: "in this project, create a ruby app
      // to deliver highest quality ingredients based on user's loc".
      args = buildArgs({ city, lat, lng, limit, ingredientsPath, citiesPath });
    } catch (error) {
      reject(error);
      return;
    }

    execFile('ruby', [binary, ...args], { timeout: timeoutMs, maxBuffer: 5 * 1024 * 1024 }, (err, stdout, stderr) => {
      if (err) {
        const details = stderr ? stderr.toString() : err.message;
        reject(new Error(`Ingredient delivery failed: ${details}`));
        return;
      }

      try {
        const parsed = JSON.parse(stdout);
        resolve(parsed);
      } catch (parseError) {
        reject(new Error(`Invalid ingredient delivery response: ${parseError.message}`));
      }
    });
  });
}

module.exports = {
  fetchHighestQualityIngredients
};
