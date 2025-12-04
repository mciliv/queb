/**
 * UI Context Automation
 * Tracks and automates UI context detection for better AI assistance
 */

class UIContextAutomation {
  constructor() {
    this.changeHistory = [];
    this.maxHistorySize = 20;
  }

  getResponsiveBreakpoint() {
    const width = window.innerWidth;
    if (width < 600) return 'mobile';
    if (width < 768) return 'tablet';
    return 'desktop';
  }

  isSignificantChange(mutations) {
    return mutations.some(mutation => {
      if (mutation.type === 'childList' && mutation.addedNodes.length > 0) {
        return true;
      }
      if (mutation.type === 'attributes') {
        const significantAttrs = ['class', 'style', 'data-state'];
        return significantAttrs.includes(mutation.attributeName);
      }
      return false;
    });
  }

  recordChange(description) {
    this.changeHistory.push({
      timestamp: new Date().toISOString(),
      description
    });
    
    // Keep only last 20 changes
    if (this.changeHistory.length > this.maxHistorySize) {
      this.changeHistory.shift();
    }
  }

  debounce(fn, delay) {
    let timeoutId;
    return function(...args) {
      clearTimeout(timeoutId);
      timeoutId = setTimeout(() => fn.apply(this, args), delay);
    };
  }

  formatContextForPrompt(context) {
    const lines = [];
    
    if (context.current_ui_state) {
      const state = context.current_ui_state;
      lines.push(`## Current UI State (Auto-Generated ${state.timestamp})`);
      lines.push('');
      
      if (state.viewport) {
        lines.push(`**Viewport**: ${state.viewport.width}x${state.viewport.height} (${state.responsiveState})`);
      }
      
      if (state.layout) {
        lines.push(`**Layout**: Sidebar ${state.layout.sidebarVisible ? 'visible' : 'hidden'}, Safe area right: ${state.layout.safeAreaRight}`);
      }
      
      if (state.activeElements) {
        lines.push(`**Active Elements**: Focus: ${state.activeElements.focusedElement}, Visible sections: ${state.activeElements.visibleSections.join(', ')}`);
      }
      
      if (state.overlappingElements && state.overlappingElements.length > 0) {
        lines.push(`**Overlapping Elements**: ${state.overlappingElements.length} found`);
      }
      
      lines.push('');
    }
    
    if (context.rule_compliance) {
      lines.push(`**Rule Compliance Score**: ${context.rule_compliance.score}%`);
      lines.push('');
    }
    
    if (context.active_issues && context.active_issues.length > 0) {
      lines.push('### Active Issues');
      context.active_issues.forEach(issue => {
        lines.push(`- [${issue.priority.toUpperCase()}] ${issue.type}: ${issue.description}`);
      });
      lines.push('');
    }
    
    if (context.suggested_focus && context.suggested_focus.length > 0) {
      lines.push('### Suggested Focus');
      context.suggested_focus.forEach(focus => {
        lines.push(`- ${focus}`);
      });
      lines.push('');
    }
    
    if (context.recent_changes && context.recent_changes.length > 0) {
      lines.push('### Recent Changes');
      context.recent_changes.forEach(change => {
        lines.push(`- ${change.timestamp}: ${change.description}`);
      });
      lines.push('');
    }
    
    return lines.join('\n');
  }
}

// Export for browser environment
if (typeof window !== 'undefined') {
  window.uiContextAutomation = new UIContextAutomation();
}

// Export for Node/test environment
if (typeof module !== 'undefined' && module.exports) {
  module.exports = UIContextAutomation;
}

