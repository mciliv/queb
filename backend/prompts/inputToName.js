// inputToName: build a prompt that converts arbitrary text input → canonical molecule names (optionally with PubChem CID)

const { buildNamesOnlyPrompt } = require('./names-only');

function inputToName(objectDescription) {
  return buildNamesOnlyPrompt(objectDescription);
}

module.exports = { inputToName };


