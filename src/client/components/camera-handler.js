/**
 * Camera Handler - Mock implementation for testing
 */

export const cameraHandler = {
  initialize: jest.fn(),
  startCamera: jest.fn(),
  stopCamera: jest.fn(),
  captureImage: jest.fn(),
  switchCamera: jest.fn(),
  getAvailableCameras: jest.fn(() => Promise.resolve([])),
  hasPermission: jest.fn(() => Promise.resolve(true)),
  requestPermission: jest.fn(() => Promise.resolve(true))
};

export default cameraHandler;
