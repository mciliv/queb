import logger from './logger.js';

export const KEYBOARD_SHORTCUTS = {
  FOCUS_INPUT: {
    key: 'k',
    modifiers: { meta: true, ctrl: true, alt: false, shift: false },
    action: 'focusInput',
    description: 'Focus input field',
    targetId: 'object-input'
  },
  CAMERA_MODE: {
    key: '1',
    modifiers: { meta: true, ctrl: true, alt: false, shift: true },
    action: 'cameraMode',
    description: 'Switch to camera mode'
  },
  PHOTO_MODE: {
    key: '2',
    modifiers: { meta: true, ctrl: true, alt: false, shift: true },
    action: 'photoMode',
    description: 'Switch to photo mode'
  },
  LINK_MODE: {
    key: '3',
    modifiers: { meta: true, ctrl: true, alt: false, shift: true },
    action: 'linkMode',
    description: 'Switch to link mode'
  }
};

export function createKeyboardHandler(actions = {}) {
  const platformIsMac = typeof navigator !== 'undefined' && (/Mac/i.test(navigator.platform) || /Mac/i.test(navigator.userAgent));

  return (keyboardEvent) => {
    logger.info(`🎹 Keyboard shortcut attempt:`, {
      key: keyboardEvent.key,
      metaKey: keyboardEvent.metaKey,
      ctrlKey: keyboardEvent.ctrlKey,
      altKey: keyboardEvent.altKey,
      shiftKey: keyboardEvent.shiftKey,
    });
    
    const deviceIsMobile = typeof window !== 'undefined' && window.matchMedia?.('(pointer: coarse)').matches;
    if (deviceIsMobile) {
      return;
    }

    const targetElement = keyboardEvent.target;
    const userIsTypingInFormField = targetElement.tagName === 'INPUT' || targetElement.tagName === 'TEXTAREA' || targetElement.isContentEditable;
    if (userIsTypingInFormField) {
      return;
    }

    for (const [shortcutName, shortcutDefinition] of Object.entries(KEYBOARD_SHORTCUTS)) {
      const { key, modifiers, targetId, action } = shortcutDefinition;

      // Check if Cmd (Mac) or Ctrl (Windows/Linux) is pressed when the shortcut requires it
      const requiresCtrlOrMeta = modifiers.meta || modifiers.ctrl;
      const ctrlOrMetaKeyPressed = requiresCtrlOrMeta ? (platformIsMac ? keyboardEvent.metaKey : keyboardEvent.ctrlKey) : true;

      const allModifiersMatch =
        ctrlOrMetaKeyPressed &&
        keyboardEvent.altKey === !!modifiers.alt &&
        keyboardEvent.shiftKey === !!modifiers.shift;

      const keyMatchesAndModifiersMatch = keyboardEvent.key.toLowerCase() === key.toLowerCase() && allModifiersMatch;
      
      if (keyMatchesAndModifiersMatch) {
        logger.debug(`🐛 Shortcut match found: ${shortcutName}`, { shortcut: shortcutDefinition });
        keyboardEvent.preventDefault();
        keyboardEvent.stopPropagation();

        const shortcutTargetsElement = !!targetId;
        if (shortcutTargetsElement) {
          const targetElementInDOM = document.getElementById(targetId);
          logger.debug(`Looking for element with ID: ${targetId}, found:`, targetElementInDOM);
          
          const elementExists = !!targetElementInDOM;
          if (elementExists) {
            targetElementInDOM.focus();
            logger.info(`✅ Shortcut ${shortcutName} executed: focused ${targetId}`);
          } else {
            logger.warn(`❌ Shortcut ${shortcutName} failed: element ${targetId} not found`);
          }
        } else {
          const actionExists = !!actions[action];
          if (actionExists) {
            actions[action]();
            logger.info(`✅ Shortcut ${shortcutName} executed: ${action}`);
          } else {
            logger.warn(`❌ Shortcut ${shortcutName} failed: action ${action} not found`);
          }
        }
        return;
      }
    }
  };
}

export function validateAllShortcutsConfigured() {
  const requiredShortcutNames = ['FOCUS_INPUT', 'CAMERA_MODE', 'PHOTO_MODE', 'LINK_MODE'];
  
  for (const shortcutName of requiredShortcutNames) {
    const shortcutExists = !!KEYBOARD_SHORTCUTS[shortcutName];
    if (!shortcutExists) {
      logger.error(`Missing required shortcut: ${shortcutName}`);
      return false;
    }
    
    const shortcutConfiguration = KEYBOARD_SHORTCUTS[shortcutName];
    const configurationIsComplete = shortcutConfiguration.key && shortcutConfiguration.action && shortcutConfiguration.modifiers;
    if (!configurationIsComplete) {
      logger.error(`Invalid shortcut configuration for ${shortcutName}`);
      return false;
    }
  }
  return true;
}