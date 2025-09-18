// UI Context Automation Module
// Provides automated UI state tracking and analysis

class UIContextAutomation {
  constructor() {
    this.changeHistory = [];
    this.maxHistorySize = 20;
    this.debouncedFunctions = new Map();
  }

  getResponsiveBreakpoint() {
    const width = window.innerWidth;
    if (width < 480) return 'mobile';
    if (width < 768) return 'tablet';
    return 'desktop';
  }

  isSignificantChange(mutations) {
    return mutations.some(mutation => {
      if (mutation.type === 'childList' && mutation.addedNodes.length > 0) return true;
      if (mutation.type === 'attributes' && mutation.attributeName === 'class') return true;
      return false;
    });
  }

  recordChange(description) {
    const change = {
      timestamp: new Date().toISOString(),
      description: description
    };

    this.changeHistory.push(change);

    // Keep only the most recent entries
    if (this.changeHistory.length > this.maxHistorySize) {
      this.changeHistory = this.changeHistory.slice(-this.maxHistorySize);
    }
  }

  debounce(func, wait) {
    const key = func.toString();
    if (this.debouncedFunctions.has(key)) {
      clearTimeout(this.debouncedFunctions.get(key));
    }

    return (...args) => {
      const timeout = setTimeout(() => func.apply(this, args), wait);
      this.debouncedFunctions.set(key, timeout);
    };
  }

  formatContextForPrompt(context) {
    let markdown = `## Current UI State (Auto-Generated ${context.current_ui_state.timestamp})\n\n`;

    markdown += `**Viewport**: ${context.current_ui_state.viewport.width}x${context.current_ui_state.viewport.height} (${context.current_ui_state.responsiveState})\n\n`;

    if (context.active_issues && context.active_issues.length > 0) {
      markdown += '### Active Issues\n';
      context.active_issues.forEach(issue => {
        markdown += `- [${issue.priority.toUpperCase()}] ${issue.type}: ${issue.description}\n`;
      });
      markdown += '\n';
    }

    if (context.suggested_focus && context.suggested_focus.length > 0) {
      markdown += '### Suggested Focus\n';
      context.suggested_focus.forEach(focus => {
        markdown += `- ${focus}\n`;
      });
      markdown += '\n';
    }

    if (context.recent_changes && context.recent_changes.length > 0) {
      markdown += '### Recent Changes\n';
      context.recent_changes.forEach(change => {
        markdown += `- ${change.description}\n`;
      });
      markdown += '\n';
    }

    return markdown;
  }
}

// Make it available globally for tests
if (typeof window !== 'undefined') {
  window.uiContextAutomation = {
    constructor: UIContextAutomation
  };
}

module.exports = UIContextAutomation;

