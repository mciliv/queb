// Manual mock for ui-utils.js
export const uiManager = {
  createColumn: jest.fn(() => ({ innerHTML: '' })),
  createLoadingColumn: jest.fn(() => ({ remove: jest.fn() })),
  createErrorMessage: jest.fn(() => ({ remove: jest.fn() })),
  fileToBase64: jest.fn().mockResolvedValue('mockbase64data'),
  urlToBase64: jest.fn().mockResolvedValue('mockbase64data')
};