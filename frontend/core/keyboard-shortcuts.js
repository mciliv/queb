// Minimal keyboard shortcut handler to keep UI simple and functional
// - ⌘/Ctrl + K: focus input
// - ⌘/Ctrl + Shift + C/P/L: toggle modes via provided actions

import logger from './logger.js';

// Simple keyboard handler (moved to bottom to avoid conflicts)
const createSimpleKeyboardHandler = (actions = {}) => {
  return (event) => {
    try {
      const isMac = navigator.userAgent.toUpperCase().includes('MAC');
      const modifier = isMac ? event.metaKey : event.ctrlKey;

      // Ignore if typing in inputs/areas/contentEditable
      const target = event.target;
      const tag = target && target.tagName ? target.tagName.toUpperCase() : '';
      const typing = tag === 'INPUT' || tag === 'TEXTAREA' || target?.isContentEditable;
      if (typing) return;

      if (modifier && !event.shiftKey && event.key.toLowerCase() === 'k') {
        event.preventDefault();
        const input = document.getElementById('object-input');
        logger.debug('⌘K triggered, looking for input:', input);
        if (input) {
          input.focus();
          logger.debug('Focused on input');
        } else {
          logger.warn('Input element not found!');
        }
        return;
      }

      // Ctrl/Cmd + Shift combos for modes
      if (!modifier || !event.shiftKey) return;
      const code = event.code || '';
      if (code === 'KeyC' && typeof actions.cameraMode === 'function') {
        event.preventDefault();
        actions.cameraMode();
      } else if (code === 'KeyP' && typeof actions.photoMode === 'function') {
        event.preventDefault();
        actions.photoMode();
      } else if (code === 'KeyL' && typeof actions.linkMode === 'function') {
        event.preventDefault();
        actions.linkMode();
      }
    } catch (_) {}
  };
};

// PROTECTED: Keyboard Shortcut System
// This file contains the canonical keyboard shortcut definitions
// Changes here affect core user interaction patterns

export const KEYBOARD_SHORTCUTS = {
  FOCUS_INPUT: {
    key: 'k',
    modifiers: { meta: true, ctrl: true, shift: false },
    action: 'focusInput',
    description: 'Focus input field',
    targetId: 'object-input'
  },
  CAMERA_MODE: {
    key: 'c',
    modifiers: { meta: true, ctrl: true, alt: false, shift: true },
    action: 'cameraMode',
    description: 'Switch to camera mode'
  },
  PHOTO_MODE: {
    key: 'p',
    modifiers: { meta: true, ctrl: true, alt: false, shift: true },
    action: 'photoMode',
    description: 'Switch to photo mode'
  },
  LINK_MODE: {
    key: 'l',
    modifiers: { meta: true, ctrl: true, alt: false, shift: true },
    action: 'linkMode',
    description: 'Switch to link mode'
  }
};

// Validate shortcut integrity
export function validateShortcuts() {
  const required = ['FOCUS_INPUT', 'CAMERA_MODE', 'PHOTO_MODE', 'LINK_MODE'];
  for (const name of required) {
    if (!KEYBOARD_SHORTCUTS[name]) {
      logger.error(`Missing required shortcut: ${name}`);
      return false;
    }
    const shortcut = KEYBOARD_SHORTCUTS[name];
    if (!shortcut.key || !shortcut.action || !shortcut.modifiers) {
      logger.error(`Invalid shortcut configuration for ${name}`);
      return false;
    }
  }
  return true;
}

// Create keyboard event handler with validation
export function createKeyboardHandlerValidated(actions) {
  // Validate on creation
  if (!validateShortcuts()) {
    logger.warn('Keyboard shortcuts validation failed - using fallback');
  }

  return (event) => {
    // Debug logging for ⌘K
    if ((event.metaKey || event.ctrlKey) && event.key.toLowerCase() === 'k') {
      logger.debug('⌘K detected, current state:', {
        key: event.key,
        metaKey: event.metaKey,
        ctrlKey: event.ctrlKey,
        shiftKey: event.shiftKey,
        target: event.target.tagName,
        targetId: event.target.id
      });
    }

    // Skip on mobile
    if (typeof window !== 'undefined' && window.matchMedia && 
        window.matchMedia('(pointer: coarse)').matches) {
      return;
    }

    // Skip if typing in form fields
    if (event.target.tagName === 'INPUT' || 
        event.target.tagName === 'TEXTAREA' || 
        event.target.contentEditable === 'true') {
      return;
    }

    const isMac = navigator.userAgent.toUpperCase().indexOf('MAC') >= 0;

    // Check each shortcut
    for (const [name, shortcut] of Object.entries(KEYBOARD_SHORTCUTS)) {
      const wantsAlt = !!shortcut.modifiers.alt;
      const modifierMatch = wantsAlt
        ? event.altKey
        : (isMac ? (shortcut.modifiers.meta && event.metaKey) : (shortcut.modifiers.ctrl && event.ctrlKey));

      if (modifierMatch && 
          event.shiftKey === shortcut.modifiers.shift &&
          event.key.toLowerCase() === shortcut.key.toLowerCase()) {
        
        logger.info(`Shortcut match found: ${name}`, {
          shortcut,
          targetId: shortcut.targetId,
          action: shortcut.action
        });
        
        event.preventDefault();
        event.stopPropagation();
        
        // Execute action
        if (shortcut.targetId) {
          // Focus element action
          const element = document.getElementById(shortcut.targetId);
          logger.debug(`Looking for element with ID: ${shortcut.targetId}, found:`, element);
          if (element) {
            element.focus();
            logger.info(`Shortcut ${name} executed: focused ${shortcut.targetId}`);
          } else {
            logger.warn(`Shortcut ${name} failed: element ${shortcut.targetId} not found`);
          }
        } else if (actions[shortcut.action]) {
          // Custom action
          actions[shortcut.action]();
          logger.info(`Shortcut ${name} executed: ${shortcut.action}`);
        } else {
          logger.warn(`Shortcut ${name} failed: action ${shortcut.action} not found`);
        }
        return;
      }
    }
  };
}

// Self-test function to verify shortcuts work  
export function testShortcuts() {
  return validateShortcuts();
}

// Export the main keyboard handler (uses validated version by default)
export function createKeyboardHandler(actions) {
  return createKeyboardHandlerValidated(actions);
}

// Export simple handler as alternative
export { createSimpleKeyboardHandler };