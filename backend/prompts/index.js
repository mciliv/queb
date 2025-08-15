// Centralized prompt exports with clear input→output aliases

const { buildChemicalAnalysisInstructions } = require('./chemical-analysis-instructions');
const { buildNamesOnlyPrompt, listChemicalsFromTextSpecification } = require('./names-only');
const { buildNameToSmilesPrompt } = require('./name-to-smiles');
const { inputToName } = require('./inputToName');
const { inputToSmiles } = require('./inputToSmiles');

// Aliases that clarify intent using input→output naming
const imageToSmiles_instructions = buildChemicalAnalysisInstructions;
const textToNames_prompt = buildNamesOnlyPrompt; // or listChemicalsFromTextSpecification
const namesToSmiles_prompt = buildNameToSmilesPrompt;

module.exports = {
  // Original exports
  buildChemicalAnalysisInstructions,
  buildNamesOnlyPrompt,
  listChemicalsFromTextSpecification,
  buildNameToSmilesPrompt,
  inputToName,
  inputToSmiles,
  // Aliases (preferred naming style)
  imageToSmiles_instructions,
  textToNames_prompt,
  namesToSmiles_prompt,
};


