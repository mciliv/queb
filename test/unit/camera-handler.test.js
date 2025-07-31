// test/unit/camera-handler.test.js - CameraHandler image analysis tests

// Mock DOM environment first
const { JSDOM } = require('jsdom');
const dom = new JSDOM('<!DOCTYPE html><html><body></body></html>');
global.window = dom.window;
global.document = dom.window.document;
global.navigator = dom.window.navigator;
global.location = dom.window.location;
global.URL = dom.window.URL;
global.CustomEvent = dom.window.CustomEvent;
global.Image = dom.window.Image;
global.atob = dom.window.atob;
global.Blob = dom.window.Blob;
global.File = dom.window.File;

// Mock fetch
global.fetch = jest.fn();

// Mock Canvas API
const mockCanvas = {
  width: 0,
  height: 0,
  getContext: jest.fn(() => ({
    drawImage: jest.fn(),
    strokeRect: jest.fn(),
    save: jest.fn(),
    restore: jest.fn(),
    imageSmoothingEnabled: false,
    strokeStyle: '',
    lineWidth: 0
  })),
  toDataURL: jest.fn(() => 'data:image/jpeg;base64,mockbase64data')
};

global.HTMLCanvasElement = jest.fn(() => mockCanvas);
document.createElement = jest.fn((tagName) => {
  if (tagName === 'canvas') return mockCanvas;
  if (tagName === 'img') {
    const img = {
      src: '',
      dataset: {},
      onload: null,
      onerror: null,
      addEventListener: jest.fn(),
      getBoundingClientRect: jest.fn(() => ({
        left: 0,
        top: 0,
        width: 300,
        height: 200
      })),
      width: 300,
      height: 200,
      naturalWidth: 300,
      naturalHeight: 200
    };
    // Simulate successful image load after a brief delay
    setTimeout(() => {
      if (img.onload) img.onload();
    }, 0);
    return img;
  }
  if (tagName === 'div') {
    return {
      className: '',
      innerHTML: '',
      textContent: '',
      style: {},
      appendChild: jest.fn(),
      addEventListener: jest.fn(),
      remove: jest.fn(),
      querySelector: jest.fn()
    };
  }
  return {
    className: '',
    innerHTML: '',
    textContent: '',
    onclick: null,
    appendChild: jest.fn(),
    addEventListener: jest.fn(),
    remove: jest.fn(),
    querySelector: jest.fn(),
    style: {}
  };
});



// Create shared mock objects that can be accessed in tests
const mockUIManager = {
  createColumn: jest.fn(() => ({ innerHTML: '' })),
  createLoadingColumn: jest.fn(() => ({ remove: jest.fn() })),
  createErrorMessage: jest.fn(() => ({ remove: jest.fn() })),
  fileToBase64: jest.fn().mockResolvedValue('mockbase64data'),
  urlToBase64: jest.fn().mockResolvedValue('mockbase64data')
};

const mockPaymentManager = {
  checkPaymentMethod: jest.fn().mockResolvedValue(true),
  incrementUsage: jest.fn().mockResolvedValue()
};

// Use __mocks__ for automatic mocking
jest.mock('../../frontend/components/ui-utils.js');
jest.mock('../../frontend/components/payment.js');

// Add additional global mocks
global.alert = jest.fn();
window.updateScrollHandles = jest.fn();
document.dispatchEvent = jest.fn();
// Create a shared object to hold element values
const mockElements = {
  'photo-url': { value: '', addEventListener: jest.fn(), remove: jest.fn() },
  'url-analyze': { addEventListener: jest.fn(), remove: jest.fn() },
  'photo-upload': { addEventListener: jest.fn(), files: [], remove: jest.fn() },
  'photo-options': { innerHTML: '', appendChild: jest.fn(), remove: jest.fn() }
};

document.getElementById = jest.fn((id) => {
  if (mockElements[id]) return mockElements[id];
  return { innerHTML: '', appendChild: jest.fn(), querySelector: jest.fn(), remove: jest.fn() };
});
document.querySelector = jest.fn(() => ({ appendChild: jest.fn() }));

// Helper to setup DOM
function setupTestDOM() {
  document.body.innerHTML = `
    <input id="photo-upload" type="file" />
    <input id="photo-url" type="text" />
    <button id="url-analyze">Analyze URL</button>
    <div id="photo-options"></div>
    <div class="snapshots-container"></div>
    <template id="photo-upload-template">
      <input id="photo-upload" type="file" />
    </template>
  `;
}

describe('CameraHandler Tests', () => {
  let CameraHandler;
  let cameraHandler;

  beforeAll(async () => {
    setupTestDOM();
    
    // Use the testable version with dependency injection
    const { CameraHandler: TestableHandler } = require('./camera-handler-testable.js');
    CameraHandler = TestableHandler;
    cameraHandler = new CameraHandler(mockUIManager, mockPaymentManager);
  });

  beforeEach(() => {
    setupTestDOM();
    jest.clearAllMocks();
    fetch.mockClear();
    
    // Reset mock elements
    mockElements['photo-url'].value = '';
    mockElements['photo-options'].innerHTML = '';
    
    // Reset mock implementations
    mockUIManager.fileToBase64.mockResolvedValue('mockbase64data');
    mockUIManager.urlToBase64.mockResolvedValue('mockbase64data');
    mockUIManager.createErrorMessage.mockReturnValue({ remove: jest.fn() });
    mockPaymentManager.checkPaymentMethod.mockResolvedValue(true);
    mockPaymentManager.incrementUsage.mockResolvedValue();
    
    // Reset fetch mock
    fetch.mockResolvedValue({
      ok: true,
      json: jest.fn().mockResolvedValue({
        chemicals: ['H2O'],
        object: 'water'
      })
    });
  });

  describe('constructor', () => {
    test('should detect mobile devices', () => {
      // Mock mobile user agent
      Object.defineProperty(navigator, 'userAgent', {
        value: 'Mozilla/5.0 (iPhone; CPU iPhone OS 14_0 like Mac OS X)',
        configurable: true
      });
      
      const mobileHandler = new CameraHandler();
      expect(mobileHandler.isMobile).toBe(true);
    });

    test('should detect desktop devices', () => {
      // Mock desktop user agent
      Object.defineProperty(navigator, 'userAgent', {
        value: 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36',
        configurable: true
      });
      
      const desktopHandler = new CameraHandler();
      expect(desktopHandler.isMobile).toBe(false);
    });
  });

  describe('handlePhotoUpload', () => {
    test('should handle valid image file upload', async () => {
      const mockFile = new File(['mock data'], 'test.jpg', { type: 'image/jpeg' });
      const event = { target: { files: [mockFile] } };
      
      const displaySpy = jest.spyOn(cameraHandler, 'displayUploadedImage').mockResolvedValue();
      
      await cameraHandler.handlePhotoUpload(event);
      
      expect(displaySpy).toHaveBeenCalledWith(mockFile);
      
      displaySpy.mockRestore();
    });

    test('should reject non-image files', async () => {
      const mockFile = new File(['mock data'], 'test.txt', { type: 'text/plain' });
      const event = { target: { files: [mockFile] } };
      
      const alertSpy = jest.spyOn(window, 'alert').mockImplementation();
      const displaySpy = jest.spyOn(cameraHandler, 'displayUploadedImage').mockResolvedValue();
      
      await cameraHandler.handlePhotoUpload(event);
      
      expect(alertSpy).toHaveBeenCalledWith('Please select an image file');
      expect(displaySpy).not.toHaveBeenCalled();
      
      alertSpy.mockRestore();
      displaySpy.mockRestore();
    });

    test('should handle empty file selection', async () => {
      const event = { target: { files: [] } };
      
      const displaySpy = jest.spyOn(cameraHandler, 'displayUploadedImage').mockResolvedValue();
      
      await cameraHandler.handlePhotoUpload(event);
      
      expect(displaySpy).not.toHaveBeenCalled();
      
      displaySpy.mockRestore();
    });
  });

  describe('handleUrlAnalysis', () => {
    test('should analyze valid image URL', async () => {
      const mockUrl = 'https://example.com/image.jpg';
      mockElements['photo-url'].value = mockUrl;
      
      const displaySpy = jest.spyOn(cameraHandler, 'displayUploadedImage').mockResolvedValue();
      
      await cameraHandler.handleUrlAnalysis();
      
      expect(mockUIManager.urlToBase64).toHaveBeenCalledWith(mockUrl);
      expect(displaySpy).toHaveBeenCalled();
      expect(mockElements['photo-url'].value).toBe('');
      
      displaySpy.mockRestore();
    });

    test('should reject empty URL', async () => {
      mockElements['photo-url'].value = '';
      
      const alertSpy = jest.spyOn(window, 'alert').mockImplementation();
      
      await cameraHandler.handleUrlAnalysis();
      
      expect(alertSpy).toHaveBeenCalledWith('Please enter an image URL');
      
      alertSpy.mockRestore();
    });

    test('should reject invalid URL format', async () => {
      // Set the mock element value directly
      mockElements['photo-url'].value = 'not-a-valid-url';
      
      const alertSpy = jest.spyOn(window, 'alert').mockImplementation();
      
      await cameraHandler.handleUrlAnalysis();
      
      expect(alertSpy).toHaveBeenCalledWith('Please enter a valid URL');
      
      alertSpy.mockRestore();
    });

    test('should handle URL loading errors', async () => {
      mockElements['photo-url'].value = 'https://example.com/nonexistent.jpg';
      
      mockUIManager.urlToBase64.mockRejectedValueOnce(new Error('Failed to load'));
      
      const errorSpy = jest.spyOn(cameraHandler, 'createClosableErrorMessage').mockImplementation();
      
      await cameraHandler.handleUrlAnalysis();
      
      expect(errorSpy).toHaveBeenCalledWith('Error loading image from URL: Failed to load');
      
      errorSpy.mockRestore();
    });
  });

  describe('displayUploadedImage', () => {
    test('should display uploaded image with mobile reticle', async () => {
      cameraHandler.isMobile = true;
      const mockFile = new File(['mock data'], 'test.jpg', { type: 'image/jpeg' });
      
      await cameraHandler.displayUploadedImage(mockFile);
      
      expect(mockUIManager.fileToBase64).toHaveBeenCalledWith(mockFile);
      expect(document.createElement).toHaveBeenCalledWith('div'); // For crosshair
      expect(document.createElement).toHaveBeenCalledWith('img');
    });

    test('should display uploaded image without mobile reticle on desktop', async () => {
      cameraHandler.isMobile = false;
      const mockFile = new File(['mock data'], 'test.jpg', { type: 'image/jpeg' });
      
      await cameraHandler.displayUploadedImage(mockFile);
      
      expect(mockUIManager.fileToBase64).toHaveBeenCalledWith(mockFile);
      expect(document.createElement).toHaveBeenCalledWith('img');
    });

    test('should handle image processing errors', async () => {
      const mockFile = new File(['mock data'], 'test.jpg', { type: 'image/jpeg' });
      mockUIManager.fileToBase64.mockRejectedValueOnce(new Error('Processing failed'));
      
      const errorSpy = jest.spyOn(cameraHandler, 'createClosableErrorMessage').mockImplementation();
      
      await cameraHandler.displayUploadedImage(mockFile);
      
      expect(errorSpy).toHaveBeenCalledWith('Error processing image: Processing failed');
      
      errorSpy.mockRestore();
    });
  });

  describe('handleImageClick', () => {
    test('should analyze image click successfully', async () => {
      const mockImg = {
        dataset: { imageBase64: 'mockbase64data' },
        getBoundingClientRect: () => ({ left: 0, top: 0, width: 300, height: 200 }),
        src: 'data:image/jpeg;base64,mockbase64data'
      };
      const mockEvent = { clientX: 150, clientY: 100 };
      
      // Mock the global Image constructor
      const originalImage = global.Image;
      global.Image = jest.fn().mockImplementation(() => {
        const img = {
          src: '',
          width: 400,
          height: 300,
          onload: null,
          onerror: null
        };
        setTimeout(() => { if (img.onload) img.onload(); }, 0);
        return img;
      });
      
      fetch.mockResolvedValueOnce({
        ok: true,
        json: jest.fn().mockResolvedValue({
          output: { object: 'Test Object', chemicals: [] }
        })
      });
      
      const emitSpy = jest.spyOn(cameraHandler, 'emitAnalysisResult').mockImplementation();
      
      await cameraHandler.handleImageClick(mockEvent, mockImg);
      
      expect(fetch).toHaveBeenCalledWith('/image-molecules', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: expect.stringContaining('"imageBase64":"mockbase64data"')
      });
      
      expect(emitSpy).toHaveBeenCalledWith(
        { output: { object: 'Test Object', chemicals: [] } },
        'Photo',
        'Unknown',
        false,
        'mockbase64data'
      );
      
      global.Image = originalImage;
      emitSpy.mockRestore();
    });

    test('should handle payment check failure', async () => {
      mockPaymentManager.checkPaymentMethod.mockRejectedValueOnce(new Error('Payment check failed'));
      
      const mockImg = { 
        dataset: { imageBase64: 'mockbase64data' },
        getBoundingClientRect: () => ({ left: 0, top: 0, width: 300, height: 200 }),
        src: 'data:image/jpeg;base64,mockbase64data'
      };
      const mockEvent = { clientX: 150, clientY: 100 };
      
      // Mock the global Image constructor
      const originalImage = global.Image;
      global.Image = jest.fn().mockImplementation(() => {
        const img = {
          src: '',
          width: 400,
          height: 300,
          onload: null,
          onerror: null
        };
        setTimeout(() => { if (img.onload) img.onload(); }, 0);
        return img;
      });
      
      fetch.mockResolvedValueOnce({
        ok: true,
        json: jest.fn().mockResolvedValue({ chemicals: [], object: 'test' })
      });
      
      await cameraHandler.handleImageClick(mockEvent, mockImg);
      
      // Should proceed despite payment error
      expect(fetch).toHaveBeenCalled();
      
      global.Image = originalImage;
    });

    test('should handle missing image data', async () => {
      const mockImg = { 
        dataset: {}, // No imageBase64 data
        getBoundingClientRect: () => ({ left: 0, top: 0, width: 300, height: 200 })
      };
      const mockEvent = { clientX: 150, clientY: 100 };
      
      const errorSpy = jest.spyOn(cameraHandler, 'createClosableErrorMessage').mockImplementation();
      
      await cameraHandler.handleImageClick(mockEvent, mockImg);
      
      expect(errorSpy).toHaveBeenCalledWith('No image data available');
      
      errorSpy.mockRestore();
    });

    test('should handle API errors gracefully', async () => {
      const mockImg = {
        dataset: { imageBase64: 'mockbase64data' },
        getBoundingClientRect: () => ({ left: 0, top: 0, width: 300, height: 200 }),
        src: 'data:image/jpeg;base64,mockbase64data'
      };
      const mockEvent = { clientX: 150, clientY: 100 };
      
      // Mock the global Image constructor
      const originalImage = global.Image;
      global.Image = jest.fn().mockImplementation(() => {
        const img = {
          src: '',
          width: 400,
          height: 300,
          onload: null,
          onerror: null
        };
        setTimeout(() => { if (img.onload) img.onload(); }, 0);
        return img;
      });
      
      fetch.mockResolvedValueOnce({
        ok: false,
        status: 500,
        statusText: 'Server Error'
      });
      
      const errorSpy = jest.spyOn(cameraHandler, 'createClosableErrorMessage').mockImplementation();
      
      await cameraHandler.handleImageClick(mockEvent, mockImg);
      
      expect(errorSpy).toHaveBeenCalledWith('Error: HTTP 500: Server Error');
      
      global.Image = originalImage;
      errorSpy.mockRestore();
    });

    test('should handle payment check errors with fallback', async () => {
      mockPaymentManager.checkPaymentMethod.mockRejectedValueOnce(new Error('Payment check failed'));
      
      const mockImg = {
        dataset: { imageBase64: 'mockbase64data' },
        getBoundingClientRect: () => ({ left: 0, top: 0, width: 300, height: 200 }),
        src: 'data:image/jpeg;base64,mockbase64data'
      };
      const mockEvent = { clientX: 150, clientY: 100 };
      
      // Mock the global Image constructor
      const originalImage = global.Image;
      global.Image = jest.fn().mockImplementation(() => {
        const img = {
          src: '',
          width: 400,
          height: 300,
          onload: null,
          onerror: null
        };
        setTimeout(() => { if (img.onload) img.onload(); }, 0);
        return img;
      });
      
      fetch.mockResolvedValueOnce({
        ok: true,
        json: jest.fn().mockResolvedValue({
          output: { object: 'Test Object', chemicals: [] }
        })
      });
      
      await cameraHandler.handleImageClick(mockEvent, mockImg);
      
      // Should proceed with analysis despite payment check error
      expect(fetch).toHaveBeenCalled();
      
      global.Image = originalImage;
    });
  });

  describe('emitAnalysisResult', () => {
    test('should emit custom event with analysis data', () => {
      const mockOutput = { object: 'Test', chemicals: [] };
      const eventSpy = jest.spyOn(document, 'dispatchEvent');
      
      cameraHandler.emitAnalysisResult(mockOutput, 'Photo', 'Test Object', false, 'base64data');
      
      expect(eventSpy).toHaveBeenCalledWith(
        expect.objectContaining({
          type: 'imageAnalysisComplete',
          detail: {
            output: mockOutput,
            icon: 'Photo',
            objectName: 'Test Object',
            useQuotes: false,
            croppedImageData: 'base64data'
          }
        })
      );
      
      eventSpy.mockRestore();
    });

    test('should emit event with default parameters', () => {
      const mockOutput = { object: 'Test', chemicals: [] };
      const eventSpy = jest.spyOn(document, 'dispatchEvent');
      
      cameraHandler.emitAnalysisResult(mockOutput, 'Photo', 'Test Object');
      
      expect(eventSpy).toHaveBeenCalledWith(
        expect.objectContaining({
          detail: expect.objectContaining({
            useQuotes: false,
            croppedImageData: null
          })
        })
      );
      
      eventSpy.mockRestore();
    });
  });

  describe('createClosableErrorMessage', () => {
    test('should create error message using UI manager', () => {
      const message = 'Test error message';
      const updateSpy = jest.spyOn(cameraHandler, 'updateScrollHandles').mockImplementation();
      
      const result = cameraHandler.createClosableErrorMessage(message);
      
      expect(mockUIManager.createErrorMessage).toHaveBeenCalledWith(
        message,
        expect.any(Object) // snapshots container
      );
      expect(updateSpy).toHaveBeenCalled();
      
      updateSpy.mockRestore();
    });
  });

  describe('updateScrollHandles', () => {
    test('should call global updateScrollHandles if available', () => {
      window.updateScrollHandles = jest.fn();
      
      cameraHandler.updateScrollHandles();
      
      expect(window.updateScrollHandles).toHaveBeenCalled();
    });

    test('should handle missing global updateScrollHandles', () => {
      delete window.updateScrollHandles;
      
      expect(() => cameraHandler.updateScrollHandles()).not.toThrow();
    });
  });

  describe('setupEventListeners', () => {
    test('should setup photo upload event listener', () => {
      const photoUpload = document.getElementById('photo-upload');
      const addEventListenerSpy = jest.spyOn(photoUpload, 'addEventListener');
      
      cameraHandler.setupEventListeners();
      
      expect(addEventListenerSpy).toHaveBeenCalledWith('change', expect.any(Function));
    });

    test('should setup URL analysis event listeners', () => {
      const photoUrl = document.getElementById('photo-url');
      const urlAnalyze = document.getElementById('url-analyze');
      
      const photoUrlSpy = jest.spyOn(photoUrl, 'addEventListener');
      const urlAnalyzeSpy = jest.spyOn(urlAnalyze, 'addEventListener');
      
      cameraHandler.setupEventListeners();
      
      expect(photoUrlSpy).toHaveBeenCalledWith('keyup', expect.any(Function));
      expect(urlAnalyzeSpy).toHaveBeenCalledWith('click', expect.any(Function));
    });

    test('should handle missing DOM elements gracefully', () => {
      // Remove elements from DOM
      document.getElementById('photo-upload').remove();
      document.getElementById('photo-url').remove();
      document.getElementById('url-analyze').remove();
      
      expect(() => cameraHandler.setupEventListeners()).not.toThrow();
    });
  });

  describe('image cropping functionality', () => {
    test('should calculate crop coordinates correctly', async () => {
      const mockImg = {
        dataset: { imageBase64: 'mockbase64data' },
        getBoundingClientRect: () => ({ left: 0, top: 0, width: 300, height: 200 }),
        src: 'data:image/jpeg;base64,mockbase64data'
      };
      const mockEvent = { clientX: 150, clientY: 100 };
      
      // Mock Image with proper onload callback
      const originalImage = global.Image;
      global.Image = jest.fn().mockImplementation(() => {
        const img = {
          onload: null,
          onerror: null,
          src: '',
          width: 600,
          height: 400
        };
        setTimeout(() => { if (img.onload) img.onload(); }, 0);
        return img;
      });
      
      fetch.mockResolvedValueOnce({
        ok: true,
        json: jest.fn().mockResolvedValue({
          output: { object: 'Test Object', chemicals: [] }
        })
      });
      
      await cameraHandler.handleImageClick(mockEvent, mockImg);
      
      // Should call fetch with calculated coordinates
      expect(fetch).toHaveBeenCalledWith('/image-molecules', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: expect.stringContaining('"x":300') // relativeX * tempImg.width
      });
      
      global.Image = originalImage;
    });
  });

  describe('error handling', () => {
    test('should handle image load errors in handleImageClick', async () => {
      const mockImg = {
        dataset: { imageBase64: 'mockbase64data' },
        getBoundingClientRect: () => ({ left: 0, top: 0, width: 300, height: 200 }),
        src: 'data:image/jpeg;base64,mockbase64data'
      };
      const mockEvent = { clientX: 150, clientY: 100 };
      
      // Mock Image with onerror callback that triggers immediately
      const originalImage = global.Image;
      global.Image = jest.fn().mockImplementation(() => {
        const img = {
          onload: null,
          onerror: null,
          src: '',
          width: 400,
          height: 300
        };
        setTimeout(() => { if (img.onerror) img.onerror(); }, 0);
        return img;
      });
      
      const errorSpy = jest.spyOn(cameraHandler, 'createClosableErrorMessage').mockImplementation();
      
      await cameraHandler.handleImageClick(mockEvent, mockImg);
      
      expect(errorSpy).toHaveBeenCalledWith('Error: Failed to process image');
      
      global.Image = originalImage;
      errorSpy.mockRestore();
    });

    test('should handle network errors in image analysis', async () => {
      const mockImg = {
        dataset: { imageBase64: 'mockbase64data' },
        getBoundingClientRect: () => ({ left: 0, top: 0, width: 300, height: 200 }),
        src: 'data:image/jpeg;base64,mockbase64data'
      };
      const mockEvent = { clientX: 150, clientY: 100 };
      
      // Mock the global Image constructor
      const originalImage = global.Image;
      global.Image = jest.fn().mockImplementation(() => {
        const img = {
          src: '',
          width: 400,
          height: 300,
          onload: null,
          onerror: null
        };
        setTimeout(() => { if (img.onload) img.onload(); }, 0);
        return img;
      });
      
      fetch.mockRejectedValueOnce(new Error('Network error'));
      
      const errorSpy = jest.spyOn(cameraHandler, 'createClosableErrorMessage').mockImplementation();
      
      await cameraHandler.handleImageClick(mockEvent, mockImg);
      
      expect(errorSpy).toHaveBeenCalledWith('Error: Network error');
      
      global.Image = originalImage;
      errorSpy.mockRestore();
    });
  });
});