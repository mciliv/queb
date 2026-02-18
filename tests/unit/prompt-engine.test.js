const promptEngine = require('../../src/core/PromptEngine');

describe('PromptEngine chemical response validation', () => {
  it('accepts chemicals with missing SMILES (null/undefined/empty) so downstream can resolve by name', () => {
    const response = {
      object: 'cacao',
      chemicals: [
        { name: 'Caffeine', smiles: null },
        { name: 'Theobromine', smiles: undefined },
        { name: 'Catechin', smiles: '' },
      ]
    };

    expect(promptEngine.validateResponse('chemical', response)).toBe(true);
  });

  it('rejects chemicals when SMILES is present but invalid', () => {
    const response = {
      object: 'cacao',
      chemicals: [
        { name: 'BadSmiles', smiles: 'NOT_A_SMILES' },
      ]
    };

    expect(promptEngine.validateResponse('chemical', response)).toBe(false);
  });

  it('rejects chemicals missing a name', () => {
    const response = {
      object: 'cacao',
      chemicals: [
        { smiles: 'O' },
      ]
    };

    expect(promptEngine.validateResponse('chemical', response)).toBe(false);
  });
});

