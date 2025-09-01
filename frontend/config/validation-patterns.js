// Validation patterns for input filtering
export const VALIDATION_PATTERNS = {
  emotions: /^(love|hate|happy|sad|angry|joy|fear|hope|dream|idea|thought|feeling|emotion)/,
  actions: /^(running|walking|talking|thinking|sleeping|eating|drinking)$/,
  shortLetters: /^[a-z]{1,2}$/,
  nonAlphabetic: /^[^a-z]*$/,
  testStrings: /^(asdf|qwerty|test|random|nothing|something|anything|everything)$/i,
  commonWords: /^(the|a|an|and|or|but|if|then|when|where|why|how|what|who)$/
};
