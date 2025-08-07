// test/unit/front-end-atom-predictor.test.js - Frontend AI integration tests
// Tests for AI analysis methods in the frontend application

// Mock DOM environment
const { JSDOM } = require('jsdom');
const dom = new JSDOM('<!DOCTYPE html><html><body></body></html>');
global.window = dom.window;
global.document = dom.window.document;
global.navigator = dom.window.navigator;
global.location = dom.window.location;

// Mock fetch globally
global.fetch = jest.fn();

// Mock 3Dmol library
global.$3Dmol = {
  createViewer: jest.fn(() => ({
    addModel: jest.fn(),
    setStyle: jest.fn(),
    zoomTo: jest.fn(),
    render: jest.fn(),
    resize: jest.fn(),
    clear: jest.fn()
  })),
  rasmolElementColors: {}
};

// Mock payment manager
const mockPaymentManager = {
  checkInitialPaymentSetup: jest.fn().mockResolvedValue(true),
  updateAccountStatus: jest.fn(),
  isDeveloperAccount: jest.fn().mockReturnValue(false),
  setupDeveloperAccount: jest.fn(),
  hidePaymentModal: jest.fn()
};

// Mock camera manager
const mockCameraManager = {
  initialize: jest.fn().mockResolvedValue(),
  isSafari: false,
  hasStoredCameraPermission: jest.fn().mockReturnValue(true),
  requestPermission: jest.fn(),
  cleanup: jest.fn()
};

// Mock camera handler
const mockCameraHandler = {
  setupEventListeners: jest.fn()
};

// Mock UI manager
const mockUIManager = {
  initialize: jest.fn(),
  setupDebuggingFunctions: jest.fn(),
  showMainApp: jest.fn(),
  switchToCameraMode: jest.fn(),
  switchToPhotoMode: jest.fn(),
  clearModeSelection: jest.fn(),
  getMoleculeName: jest.fn(chemical => chemical.name || 'Unknown'),
  createErrorMessage: jest.fn(() => document.createElement('div')),
  cleanup: jest.fn()
};

// Mock modules
jest.doMock('../../frontend/components/payment.js', () => ({
  paymentManager: mockPaymentManager
}));

jest.doMock('../../frontend/components/camera.js', () => ({
  cameraManager: mockCameraManager
}));

jest.doMock('../../frontend/components/camera-handler.js', () => ({
  cameraHandler: mockCameraHandler
}));

jest.doMock('../../frontend/components/ui-utils.js', () => ({
  uiManager: mockUIManager
}));

// Mock localStorage
const mockLocalStorage = {
  getItem: jest.fn(),
  setItem: jest.fn(),
  clear: jest.fn()
};
global.localStorage = mockLocalStorage;

// Helper to create DOM elements for tests
function setupTestDOM() {
  document.body.innerHTML = `
    <div class="snapshots-container"></div>
    <input id="object-input" type="text" />
    <div id="payment-modal" class="hidden"></div>
    <div id="processing-indicator" style="display: none;"></div>
    <div id="video-feed"></div>
    <input id="photo-upload" type="file" />
    <input id="photo-url" type="text" />
    <button id="url-analyze"></button>
    <div id="modal-backdrop"></div>
    <button id="payment-close-btn"></button>
    <button id="start-analyzing-btn"></button>
  `;
}

describe('Frontend AI Integration Tests', () => {
  let MolecularApp;
  let app;

  beforeAll(async () => {
    // Set up DOM
    setupTestDOM();
    
    // Import MolecularApp after mocking
    const appModule = await import('../../frontend/core/app.js');
    MolecularApp = appModule.MolecularApp;
  });

  beforeEach(() => {
    setupTestDOM();
    app = new MolecularApp();
    
    // Reset mocks
    fetch.mockClear();
    mockPaymentManager.checkInitialPaymentSetup.mockClear();
    mockPaymentManager.updateAccountStatus.mockClear();
    mockLocalStorage.getItem.mockClear();
    mockLocalStorage.setItem.mockClear();
  });

  describe('MolecularApp initialization', () => {
    test('should initialize with default properties', () => {
      expect(app.snapshots).toBe(null);
      expect(app.objectInput).toBe(null);
      expect(app.viewers).toEqual([]);
      expect(app.isProcessing).toBe(false);
      expect(app.hasPaymentSetup).toBe(false);
      expect(app.currentAnalysisType).toBe(null);
      expect(app.lastAnalysis).toBe(null);
    });

    test('should initialize components properly', async () => {
      await app.initialize();
      
      expect(app.snapshots).toBeTruthy();
      expect(app.objectInput).toBeTruthy();
      expect(mockUIManager.initialize).toHaveBeenCalled();
      expect(mockUIManager.setupDebuggingFunctions).toHaveBeenCalled();
      expect(mockUIManager.showMainApp).toHaveBeenCalled();
      expect(mockPaymentManager.checkPaymentRequired).toHaveBeenCalled();
    });

    test('should handle localhost developer mode setup', async () => {
      // Mock localhost environment
      Object.defineProperty(window, 'location', {
        value: { hostname: 'localhost' },
        writable: true
      });
      
      mockLocalStorage.getItem.mockReturnValue(null);
      mockPaymentManager.isDeveloperAccount.mockReturnValue(false);
      
      await app.initialize();
      
      expect(mockPaymentManager.checkPaymentRequired).toHaveBeenCalled();
      expect(app.hasPaymentSetup).toBe(true);
    });
  });

  describe('handleTextAnalysis', () => {
    beforeEach(() => {
      app.objectInput = document.getElementById('object-input');
      app.objectInput.value = 'ethanol';
      app.hasPaymentSetup = true;
    });

    test('should analyze text successfully', async () => {
      const mockResponse = {
        object: 'ethanol',
        chemicals: [{ name: 'Ethanol', smiles: 'CCO' }]
      };

      fetch.mockResolvedValueOnce({
        ok: true,
        json: jest.fn().mockResolvedValue(mockResponse)
      });

      const processAnalysisResultSpy = jest.spyOn(app, 'processAnalysisResult').mockImplementation();

      await app.handleTextAnalysis();

      expect(fetch).toHaveBeenCalledWith('/analyze-text', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ text: 'ethanol' }),
        signal: expect.any(AbortSignal)
      });

      expect(app.lastAnalysis).toEqual(mockResponse);
      expect(processAnalysisResultSpy).toHaveBeenCalledWith(mockResponse, null, 'ethanol', false, null);
      expect(app.objectInput.value).toBe('');
      
      processAnalysisResultSpy.mockRestore();
    });

    test('should handle API errors gracefully', async () => {
      fetch.mockResolvedValueOnce({
        ok: false,
        status: 500,
        text: jest.fn().mockResolvedValue('Internal Server Error')
      });

      const handleErrorSpy = jest.spyOn(app, 'handleError').mockImplementation();

      await app.handleTextAnalysis();

      expect(handleErrorSpy).toHaveBeenCalled();
      expect(app.isProcessing).toBe(false);
      
      handleErrorSpy.mockRestore();
    });

    test('should prevent analysis without payment setup', async () => {
      app.hasPaymentSetup = false;
      
      await app.handleTextAnalysis();
      
      expect(fetch).not.toHaveBeenCalled();
    });

    test('should prevent duplicate analysis when processing', async () => {
      app.isProcessing = true;
      
      await app.handleTextAnalysis();
      
      expect(fetch).not.toHaveBeenCalled();
    });

    test('should handle empty input gracefully', async () => {
      app.objectInput.value = '';
      
      await app.handleTextAnalysis();
      
      expect(fetch).not.toHaveBeenCalled();
    });
  });

  describe('processAnalysisResult', () => {
    beforeEach(() => {
      app.snapshots = document.querySelector('.snapshots-container');
    });

    test('should process analysis result with chemicals', async () => {
      const mockOutput = {
        object: 'water',
        chemicals: [{ name: 'Water', smiles: 'O' }]
      };

      const generateSDFsSpy = jest.spyOn(app, 'generateSDFs').mockResolvedValue();
      
      await app.processAnalysisResult(mockOutput, null, 'water', false, null);
      
      expect(generateSDFsSpy).toHaveBeenCalled();
      
      generateSDFsSpy.mockRestore();
    });

    test('should handle quoted object names', async () => {
      const mockOutput = {
        object: 'quoted object',
        chemicals: []
      };

      const createObjectColumnSpy = jest.spyOn(app, 'createObjectColumn').mockImplementation();
      
      await app.processAnalysisResult(mockOutput, null, 'test', true, null);
      
      expect(createObjectColumnSpy).toHaveBeenCalledWith(
        '"test"',
        [],
        [],
        "No displayable molecular structures found"
      );
      
      createObjectColumnSpy.mockRestore();
    });
  });

  describe('generateSDFs', () => {
    test('should generate SDF files successfully', async () => {
      const mockResponse = {
        sdfPaths: ['/sdf_files/CCO.sdf'],
        errors: [],
        skipped: []
      };

      fetch.mockResolvedValueOnce({
        ok: true,
        json: jest.fn().mockResolvedValue(mockResponse)
      });

      const createObjectColumnSpy = jest.spyOn(app, 'createObjectColumn').mockImplementation();

      await app.generateSDFs(['CCO'], 'ethanol', null, [{ name: 'Ethanol', smiles: 'CCO' }], null);

      expect(fetch).toHaveBeenCalledWith('/generate-sdfs', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ smiles: ['CCO'], overwrite: true })
      });

      expect(createObjectColumnSpy).toHaveBeenCalledWith(
        'ethanol',
        null,
        ['CCO'],
        ['/sdf_files/CCO.sdf'],
        [{ name: 'Ethanol', smiles: 'CCO' }],
        null
      );
      
      createObjectColumnSpy.mockRestore();
    });

    test('should handle SDF generation errors', async () => {
      const mockResponse = {
        sdfPaths: [],
        errors: ['CCO - Generation failed'],
        skipped: []
      };

      fetch.mockResolvedValueOnce({
        ok: true,
        json: jest.fn().mockResolvedValue(mockResponse)
      });

      const handleErrorSpy = jest.spyOn(app, 'handleError').mockImplementation();

      await app.generateSDFs(['CCO'], 'ethanol', null, [{ name: 'Ethanol', smiles: 'CCO' }], null);

      expect(handleErrorSpy).toHaveBeenCalled();
      
      handleErrorSpy.mockRestore();
    });

    test('should handle network errors in SDF generation', async () => {
      fetch.mockRejectedValueOnce(new Error('Network error'));

      const handleErrorSpy = jest.spyOn(app, 'handleError').mockImplementation();

      await app.generateSDFs(['CCO'], 'ethanol', null, [{ name: 'Ethanol', smiles: 'CCO' }], null);

      expect(handleErrorSpy).toHaveBeenCalledWith(expect.any(Error));
      
      handleErrorSpy.mockRestore();
    });
  });

  describe('render', () => {
    test('should render 3D molecule successfully', async () => {
      const mockSdfData = 'mock sdf content';
      const mockViewer = {
        addModel: jest.fn(),
        setStyle: jest.fn(),
        zoomTo: jest.fn(),
        render: jest.fn()
      };

      fetch.mockResolvedValueOnce({
        ok: true,
        text: jest.fn().mockResolvedValue(mockSdfData)
      });

      global.$3Dmol.createViewer.mockReturnValue(mockViewer);

      const container = document.createElement('div');
      const result = await app.render('/sdf_files/test.sdf', container);

      expect(fetch).toHaveBeenCalledWith('/sdf_files/test.sdf');
      expect($3Dmol.createViewer).toHaveBeenCalledWith(container, {
        defaultcolors: $3Dmol.rasmolElementColors
      });
      expect(mockViewer.addModel).toHaveBeenCalledWith(mockSdfData, 'sdf');
      expect(mockViewer.setStyle).toHaveBeenCalledWith({}, { sphere: { scale: 0.8 } });
      expect(mockViewer.zoomTo).toHaveBeenCalled();
      expect(mockViewer.render).toHaveBeenCalled();
      expect(result).toBe(mockViewer);
    });

    test('should handle render errors gracefully', async () => {
      fetch.mockRejectedValueOnce(new Error('Failed to load SDF'));

      const container = document.createElement('div');
      const result = await app.render('/sdf_files/invalid.sdf', container);

      expect(container.textContent).toContain('Error:');
      expect(result).toBe(null);
    });

    test('should handle HTTP errors in render', async () => {
      fetch.mockResolvedValueOnce({
        ok: false,
        status: 404
      });

      const container = document.createElement('div');
      const result = await app.render('/sdf_files/notfound.sdf', container);

      expect(container.textContent).toContain('Error:');
      expect(result).toBe(null);
    });
  });

  describe('UI state management', () => {
    test('should show processing indicator', () => {
      const indicator = document.getElementById('processing-indicator');
      
      app.showProcessing();
      
      expect(indicator.style.display).toBe('block');
    });

    test('should hide processing indicator', () => {
      const indicator = document.getElementById('processing-indicator');
      indicator.style.display = 'block';
      
      app.hideProcessing();
      
      expect(indicator.style.display).toBe('none');
    });

    test('should handle error messages', () => {
      const error = new Error('Test error message');
      
      const createClosableErrorMessageSpy = jest.spyOn(app, 'createClosableErrorMessage').mockImplementation();
      
      app.handleError(error);
      
      expect(createClosableErrorMessageSpy).toHaveBeenCalledWith('Test error message');
      
      createClosableErrorMessageSpy.mockRestore();
    });

    test('should create error messages', () => {
      app.snapshots = document.querySelector('.snapshots-container');
      
      const errorDiv = app.createClosableErrorMessage('Test error');
      
      expect(mockUIManager.createErrorMessage).toHaveBeenCalledWith('Test error', app.snapshots);
    });
  });

  describe('cleanup', () => {
    test('should cleanup all resources', () => {
      const mockViewer = { clear: jest.fn() };
      app.viewers = [mockViewer];
      
      app.cleanup();
      
      expect(mockCameraManager.cleanup).toHaveBeenCalled();
      expect(mockUIManager.cleanup).toHaveBeenCalled();
      expect(mockViewer.clear).toHaveBeenCalled();
      expect(app.viewers).toEqual([]);
    });

    test('should handle viewers without clear method', () => {
      const mockViewer = {}; // No clear method
      app.viewers = [mockViewer];
      
      expect(() => app.cleanup()).not.toThrow();
      expect(app.viewers).toEqual([]);
    });
  });

  describe('event listeners', () => {
    test('should setup text analysis event listeners', () => {
      const textInput = document.getElementById('object-input');
      const enterEvent = new KeyboardEvent('keydown', { key: 'Enter' });
      
      const handleTextAnalysisSpy = jest.spyOn(app, 'handleTextAnalysis').mockImplementation();
      
      app.setupTextAnalysis();
      textInput.dispatchEvent(enterEvent);
      
      expect(handleTextAnalysisSpy).toHaveBeenCalled();
      
      handleTextAnalysisSpy.mockRestore();
    });

    test('should setup keyboard shortcuts', async () => {
      await app.initialize();
      
      const textInput = document.getElementById('object-input');
      const focusSpy = jest.spyOn(textInput, 'focus');
      const selectSpy = jest.spyOn(textInput, 'select');
      
      // Simulate Cmd+K
      const cmdKEvent = new KeyboardEvent('keydown', { 
        key: 'k', 
        metaKey: true,
        bubbles: true 
      });
      
      document.dispatchEvent(cmdKEvent);
      
      expect(focusSpy).toHaveBeenCalled();
      expect(selectSpy).toHaveBeenCalled();
      
      focusSpy.mockRestore();
      selectSpy.mockRestore();
    });
  });
});