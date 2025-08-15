// inputToSmiles: build a prompt that converts a list of molecule names/CIDs → isomeric SMILES

const { buildNameToSmilesPrompt } = require('./name-to-smiles');

function inputToSmiles(payload) {
  return buildNameToSmilesPrompt(payload);
}

module.exports = { inputToSmiles };


