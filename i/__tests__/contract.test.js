const {
  loadOpenAIResponseContract,
  createSuccessfulAIResponse
} = require('../structuralize-api-contract');

describe('OpenAI response contract helpers', () => {
  it('loads contract from openai_respone.json', () => {
    const contract = loadOpenAIResponseContract();
    expect(contract).toEqual(
      expect.objectContaining({
        object: 'response'
      })
    );
  });

  it('creates a mock response with output_text', () => {
    const mock = createSuccessfulAIResponse('coffee', [
      { name: 'Caffeine', smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' }
    ]);

    expect(mock).toEqual(
      expect.objectContaining({
        object: 'response',
        output: [
          expect.objectContaining({
            type: 'message',
            content: [
              expect.objectContaining({
                type: 'output_text'
              })
            ]
          })
        ]
      })
    );
  });
});
