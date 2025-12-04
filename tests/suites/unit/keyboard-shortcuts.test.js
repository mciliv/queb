// Test file for keyboard shortcuts functionality
// Tests that keyboard shortcuts perform their intended actions

/**
 * @jest-environment jsdom
 */

// Mock logger to avoid console spam during tests
jest.mock('../../../src/client/logger.js', () => ({
  debug: jest.fn(),
  info: jest.fn(),
  warn: jest.fn(),
  error: jest.fn()
}));

// Import the keyboard shortcuts module
const {
  KEYBOARD_SHORTCUTS,
  createKeyboardHandler,
  validateAllShortcutsConfigured
} = require('../../../src/client/keyboard-shortcuts.js');

describe('Keyboard Shortcuts Functionality Tests', () => {
  let mockActions;
  let keyboardHandler;

  beforeEach(() => {
    // Setup DOM
    document.body.innerHTML = `
      <div>
        <input type="text" id="object-input" />
        <input type="text" id="other-input" />
        <textarea id="text-area"></textarea>
        <div contenteditable="true" id="editable-div"></div>
      </div>
    `;

    // Create mock actions
    mockActions = {
      focusInput: jest.fn(),
      cameraMode: jest.fn(),
      photoMode: jest.fn(),
      linkMode: jest.fn()
    };

    // Reset navigator.userAgent to Mac for consistent testing
    Object.defineProperty(navigator, 'userAgent', {
      writable: true,
      value: 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)'
    });

    // Create keyboard handler after setting userAgent
    keyboardHandler = createKeyboardHandler(mockActions);

    // Mock window.matchMedia to always return non-mobile (desktop)
    Object.defineProperty(window, 'matchMedia', {
      writable: true,
      value: jest.fn().mockImplementation(query => ({
        matches: false, // Always return false (desktop)
        media: query,
        onchange: null,
        addListener: jest.fn(),
        removeListener: jest.fn(),
        addEventListener: jest.fn(),
        removeEventListener: jest.fn(),
        dispatchEvent: jest.fn(),
      })),
    });
  });

  afterEach(() => {
    jest.clearAllMocks();
    
    // Reset any mocked window properties
    if (window.matchMedia && window.matchMedia.mockRestore) {
      window.matchMedia.mockRestore();
    }
    
    // Clean up DOM
    document.body.innerHTML = '';
  });

  describe('Shortcut Configuration', () => {
    test('All required shortcuts are defined', () => {
      expect(KEYBOARD_SHORTCUTS.FOCUS_INPUT).toBeDefined();
      expect(KEYBOARD_SHORTCUTS.CAMERA_MODE).toBeDefined();
      expect(KEYBOARD_SHORTCUTS.PHOTO_MODE).toBeDefined();
      expect(KEYBOARD_SHORTCUTS.LINK_MODE).toBeDefined();
    });

    test('Shortcuts have correct structure', () => {
      Object.values(KEYBOARD_SHORTCUTS).forEach(shortcut => {
        expect(shortcut).toHaveProperty('key');
        expect(shortcut).toHaveProperty('modifiers');
        expect(shortcut).toHaveProperty('action');
        expect(shortcut).toHaveProperty('description');
        expect(shortcut.modifiers).toHaveProperty('meta');
        expect(shortcut.modifiers).toHaveProperty('ctrl');
        expect(shortcut.modifiers).toHaveProperty('shift');
      });
    });

    test('Validation function works correctly', () => {
      expect(validateAllShortcutsConfigured()).toBe(true);
    });
  });

  describe('Focus Input Shortcut (Cmd/Ctrl+K)', () => {
    test('focuses input field when Cmd+K is pressed on Mac', () => {
      const input = document.getElementById('object-input');
      const focusSpy = jest.spyOn(input, 'focus');

      const event = new KeyboardEvent('keydown', {
        key: 'k',
        metaKey: true,
        ctrlKey: false,
        shiftKey: false,
        bubbles: true
      });
      
      // Set target to document body (not an input field)
      Object.defineProperty(event, 'target', { value: document.body, writable: false });

      keyboardHandler(event);
      expect(focusSpy).toHaveBeenCalled();
    });

    test('focuses input field when Ctrl+K is pressed on Windows', () => {
      // Mock Windows user agent
      Object.defineProperty(navigator, 'userAgent', {
        writable: true,
        value: 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'
      });

      // Recreate handler with Windows userAgent
      const windowsKeyboardHandler = createKeyboardHandler(mockActions);

      const input = document.getElementById('object-input');
      const focusSpy = jest.spyOn(input, 'focus');

      const event = new KeyboardEvent('keydown', {
        key: 'k',
        metaKey: false,
        ctrlKey: true,
        shiftKey: false,
        bubbles: true
      });

      // Set target to document body (not an input field)
      Object.defineProperty(event, 'target', { value: document.body, writable: false });

      windowsKeyboardHandler(event);
      expect(focusSpy).toHaveBeenCalled();
    });

    test('does not trigger when Shift is also pressed', () => {
      const input = document.getElementById('object-input');
      const focusSpy = jest.spyOn(input, 'focus');

      const event = new KeyboardEvent('keydown', {
        key: 'k',
        metaKey: true,
        shiftKey: true,
        bubbles: true
      });
      
      // Set target to document body (not an input field)
      Object.defineProperty(event, 'target', { value: document.body, writable: false });

      keyboardHandler(event);
      expect(focusSpy).not.toHaveBeenCalled();
    });
  });

  describe('Mode Toggle Shortcuts (Cmd/Ctrl+Shift+1/2/3)', () => {
    test('camera mode shortcut (Cmd+Shift+1) triggers cameraMode action', () => {
      const event = new KeyboardEvent('keydown', {
        key: '1',
        metaKey: true,
        ctrlKey: false,
        altKey: false,
        shiftKey: true,
        bubbles: true
      });
      
      // Set target to document body (not an input field)
      Object.defineProperty(event, 'target', { value: document.body, writable: false });

      keyboardHandler(event);
      expect(mockActions.cameraMode).toHaveBeenCalled();
    });

    test('photo mode shortcut (Cmd+Shift+2) triggers photoMode action', () => {
      const event = new KeyboardEvent('keydown', {
        key: '2',
        metaKey: true,
        ctrlKey: false,
        altKey: false,
        shiftKey: true,
        bubbles: true
      });
      
      // Set target to document body (not an input field)
      Object.defineProperty(event, 'target', { value: document.body, writable: false });

      keyboardHandler(event);
      expect(mockActions.photoMode).toHaveBeenCalled();
    });

    test('link mode shortcut (Cmd+Shift+3) triggers linkMode action', () => {
      const event = new KeyboardEvent('keydown', {
        key: '3',
        metaKey: true,
        ctrlKey: false,
        altKey: false,
        shiftKey: true,
        bubbles: true
      });
      
      // Set target to document body (not an input field)
      Object.defineProperty(event, 'target', { value: document.body, writable: false });

      keyboardHandler(event);
      expect(mockActions.linkMode).toHaveBeenCalled();
    });

    test('mode shortcuts work with Ctrl on Windows', () => {
      // Mock Windows user agent
      Object.defineProperty(navigator, 'userAgent', {
        writable: true,
        value: 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'
      });

      // Recreate handler with Windows userAgent
      const windowsKeyboardHandler = createKeyboardHandler(mockActions);

      const event = new KeyboardEvent('keydown', {
        key: '1',
        ctrlKey: true,
        altKey: false,
        shiftKey: true,
        bubbles: true
      });

      // Set target to document body (not an input field)
      Object.defineProperty(event, 'target', { value: document.body, writable: false });

      windowsKeyboardHandler(event);
      expect(mockActions.cameraMode).toHaveBeenCalled();
    });
  });

  describe('Shortcut Blocking in Input Fields', () => {
    test('shortcuts are ignored when typing in input field', () => {
      const input = document.getElementById('object-input');
      input.focus();

      const event = new KeyboardEvent('keydown', {
        key: 'k',
        metaKey: true,
        target: input,
        bubbles: true
      });

      // Simulate the event targeting the input
      Object.defineProperty(event, 'target', { value: input });

      keyboardHandler(event);
      expect(mockActions.focusInput).not.toHaveBeenCalled();
    });

    test('shortcuts are ignored when typing in textarea', () => {
      const textarea = document.getElementById('text-area');
      textarea.focus();

      const event = new KeyboardEvent('keydown', {
        key: 'c',
        code: 'KeyC',
        metaKey: true,
        shiftKey: true,
        target: textarea,
        bubbles: true
      });

      Object.defineProperty(event, 'target', { value: textarea });

      keyboardHandler(event);
      expect(mockActions.cameraMode).not.toHaveBeenCalled();
    });

    test('shortcuts are ignored in contentEditable elements', () => {
      const editableDiv = document.getElementById('editable-div');
      editableDiv.focus();
      
      // Ensure contentEditable is actually set to 'true' string
      editableDiv.contentEditable = 'true';

      const event = new KeyboardEvent('keydown', {
        key: 'p',
        code: 'KeyP',
        metaKey: true,
        shiftKey: true,
        bubbles: true
      });

      // Set target to the contentEditable div
      Object.defineProperty(event, 'target', { value: editableDiv, writable: false });

      keyboardHandler(event);
      expect(mockActions.photoMode).not.toHaveBeenCalled();
    });
  });

  describe('Event Prevention and Propagation', () => {
    test('prevents default behavior for valid shortcuts', () => {
      const event = new KeyboardEvent('keydown', {
        key: 'k',
        metaKey: true,
        bubbles: true
      });
      
      // Set target to document body (not an input field)
      Object.defineProperty(event, 'target', { value: document.body, writable: false });
      
      const preventDefaultSpy = jest.spyOn(event, 'preventDefault');
      
      keyboardHandler(event);
      expect(preventDefaultSpy).toHaveBeenCalled();
    });

    test('stops propagation for valid shortcuts', () => {
      const event = new KeyboardEvent('keydown', {
        key: '1',
        metaKey: true,
        ctrlKey: false,
        altKey: false,
        shiftKey: true,
        bubbles: true
      });
      
      // Set target to document body (not an input field)
      Object.defineProperty(event, 'target', { value: document.body, writable: false });
      
      const stopPropagationSpy = jest.spyOn(event, 'stopPropagation');
      
      keyboardHandler(event);
      expect(stopPropagationSpy).toHaveBeenCalled();
    });
  });

  describe('Edge Cases and Error Handling', () => {
    test('handles missing input element gracefully', () => {
      // Remove the input element
      document.getElementById('object-input').remove();

      const event = new KeyboardEvent('keydown', {
        key: 'k',
        metaKey: true,
        bubbles: true
      });
      
      // Set target to document body (not an input field)
      Object.defineProperty(event, 'target', { value: document.body, writable: false });

      expect(() => keyboardHandler(event)).not.toThrow();
    });

    test('handles missing action functions gracefully', () => {
      const handlerWithoutActions = createKeyboardHandler({});

      const event = new KeyboardEvent('keydown', {
        key: 'c',
        code: 'KeyC',
        metaKey: true,
        shiftKey: true,
        bubbles: true
      });
      
      // Set target to document body (not an input field)
      Object.defineProperty(event, 'target', { value: document.body, writable: false });

      expect(() => handlerWithoutActions(event)).not.toThrow();
    });

    test('handles events without code property', () => {
      const event = new KeyboardEvent('keydown', {
        key: 'c',
        metaKey: true,
        shiftKey: true,
        bubbles: true
      });
      
      // Set target to document body (not an input field)
      Object.defineProperty(event, 'target', { value: document.body, writable: false });
      
      // Remove code property to simulate older browsers
      delete event.code;

      expect(() => keyboardHandler(event)).not.toThrow();
    });

    test('ignores shortcuts on mobile devices', () => {
      // Override the default mock to simulate mobile device
      window.matchMedia.mockImplementation(query => ({
        matches: query === '(pointer: coarse)', // Return true for mobile detection
        media: query,
        onchange: null,
        addListener: jest.fn(),
        removeListener: jest.fn(),
        addEventListener: jest.fn(),
        removeEventListener: jest.fn(),
        dispatchEvent: jest.fn(),
      }));

      const event = new KeyboardEvent('keydown', {
        key: 'k',
        metaKey: true,
        bubbles: true
      });
      
      // Set target to document body (not an input field)
      Object.defineProperty(event, 'target', { value: document.body, writable: false });

      const input = document.getElementById('object-input');
      const focusSpy = jest.spyOn(input, 'focus');

      keyboardHandler(event);
      expect(focusSpy).not.toHaveBeenCalled();
    });
  });

  describe('Case Insensitive Key Matching', () => {
    test('shortcuts work with uppercase keys', () => {
      const input = document.getElementById('object-input');
      expect(input).toBeTruthy(); // Ensure input exists
      const focusSpy = jest.spyOn(input, 'focus');

      const event = new KeyboardEvent('keydown', {
        key: 'K', // uppercase
        metaKey: true,
        ctrlKey: false,
        shiftKey: false,
        bubbles: true
      });
      
      // Set target to document body (not an input field)
      Object.defineProperty(event, 'target', { value: document.body, writable: false });

      keyboardHandler(event);
      expect(focusSpy).toHaveBeenCalled();
    });

    test('number shortcuts work with Shift modifier', () => {
      const event = new KeyboardEvent('keydown', {
        key: '2',
        metaKey: true,
        ctrlKey: false,
        altKey: false,
        shiftKey: true,
        bubbles: true
      });
      
      // Set target to document body (not an input field)
      Object.defineProperty(event, 'target', { value: document.body, writable: false });

      keyboardHandler(event);
      expect(mockActions.photoMode).toHaveBeenCalled();
    });
  });
});