/**
 * UI Utils - Mock implementation for testing
 */

export const uiManager = {
  showLoading: jest.fn(),
  hideLoading: jest.fn(),
  showError: jest.fn(),
  showSuccess: jest.fn(),
  updateUI: jest.fn(),
  resetUI: jest.fn()
};

export default uiManager;
