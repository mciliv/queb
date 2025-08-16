// Centralized prompt exports with clear input→output aliases

const { buildChemicalAnalysisInstructions } = require('./inputToSmiles');
const { buildNamesOnlyPrompt, listChemicalsFromTextSpecification } = require('./names-only');
const { inputToName } = require('./inputToName');

// Aliases that clarify intent using input→output naming
const imageToSmiles_instructions = buildChemicalAnalysisInstructions;
const textToNames_prompt = buildNamesOnlyPrompt; // or listChemicalsFromTextSpecification
// Preferred generic grouping: inputTo*
const inputToSmilesFromImage = buildChemicalAnalysisInstructions;
const inputToNames = inputToName;
// LLM name→SMILES prompt deprecated in favor of programmatic resolution

module.exports = {
  // Original exports
  buildChemicalAnalysisInstructions,
  buildNamesOnlyPrompt,
  listChemicalsFromTextSpecification,
  inputToName,
  // inputToSmiles intentionally omitted
  // Aliases (preferred naming style)
  imageToSmiles_instructions,
  textToNames_prompt,
  inputToSmilesFromImage,
  inputToNames,
  // namesToSmiles_prompt intentionally omitted
};


