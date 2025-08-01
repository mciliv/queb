class UIContextAutomation {
  constructor() {
    this.snapshots = new Map();
    this.changeHistory = [];
    this.enabled = false;
  }

  enable() {
    this.enabled = true;
    this.setupMutationObserver();
    this.setupEventListeners();
    
  }

  generateUISnapshot() {
    return {
      timestamp: new Date().toISOString(),
      viewport: {
        width: window.innerWidth,
        height: window.innerHeight,
        devicePixelRatio: window.devicePixelRatio
      },
      layout: {
        sidebarVisible: !document.getElementById('payment-section')?.classList.contains('hidden'),
        accountButtonPosition: this.getElementPosition('.account-link'),
        textInputMargin: getComputedStyle(document.getElementById('object-input'))?.marginRight,
        safeAreaRight: getComputedStyle(document.documentElement).getPropertyValue('--safe-area-right')
      },
      activeElements: {
        focusedElement: document.activeElement?.id || null,
        visibleSections: this.getVisibleSections(),
        errorStates: this.getErrorStates()
      },
      cssVariables: this.getCSSVariables(),
      responsiveState: this.getResponsiveBreakpoint(),
      overlappingElements: this.detectOverlaps()
    };
  }

  getElementPosition(selector) {
    const element = document.querySelector(selector);
    if (!element) return null;
    const rect = element.getBoundingClientRect();
    return {
      top: rect.top,
      left: rect.left,
      width: rect.width,
      height: rect.height
    };
  }

  getVisibleSections() {
    const sections = ['camera-container', 'photo-options', 'payment-section', 'results-section'];
    return sections.filter(id => {
      const element = document.getElementById(id);
      return element && !element.classList.contains('hidden') && 
             getComputedStyle(element).display !== 'none';
    });
  }

  getErrorStates() {
    const errorElements = document.querySelectorAll('.error-message, .error-text');
    return Array.from(errorElements)
      .filter(el => el.textContent.trim())
      .map(el => ({
        element: el.className,
        message: el.textContent.trim()
      }));
  }

  getCSSVariables() {
    const computedStyle = getComputedStyle(document.documentElement);
    const variables = {};
    
    // Get all our custom properties
    const customProps = [
      '--space-xs', '--space-sm', '--space-md', '--space-lg', '--space-xl', '--space-xxl',
      '--sidebar-width', '--sidebar-gap', '--safe-area-right', '--safe-area-top',
      '--z-base', '--z-content', '--z-sidebar', '--z-modal', '--z-tooltip'
    ];
    
    customProps.forEach(prop => {
      variables[prop] = computedStyle.getPropertyValue(prop).trim();
    });
    
    return variables;
  }

  getResponsiveBreakpoint() {
    const width = window.innerWidth;
    if (width <= 480) return 'mobile';
    if (width <= 768) return 'tablet';
    return 'desktop';
  }

  detectOverlaps() {
    const overlaps = [];
    const elements = document.querySelectorAll('.account-link, #object-input, .sidebar-container, .analysis-section');
    
    for (let i = 0; i < elements.length; i++) {
      for (let j = i + 1; j < elements.length; j++) {
        if (this.elementsOverlap(elements[i], elements[j])) {
          overlaps.push({
            element1: elements[i].className || elements[i].id,
            element2: elements[j].className || elements[j].id,
            severity: this.calculateOverlapSeverity(elements[i], elements[j])
          });
        }
      }
    }
    
    return overlaps;
  }

  elementsOverlap(el1, el2) {
    const rect1 = el1.getBoundingClientRect();
    const rect2 = el2.getBoundingClientRect();
    
    return !(rect1.right < rect2.left || 
             rect1.left > rect2.right || 
             rect1.bottom < rect2.top || 
             rect1.top > rect2.bottom);
  }

  calculateOverlapSeverity(el1, el2) {
    const rect1 = el1.getBoundingClientRect();
    const rect2 = el2.getBoundingClientRect();
    
    const overlapArea = Math.max(0, Math.min(rect1.right, rect2.right) - Math.max(rect1.left, rect2.left)) *
                       Math.max(0, Math.min(rect1.bottom, rect2.bottom) - Math.max(rect1.top, rect2.top));
    
    const totalArea = (rect1.width * rect1.height) + (rect2.width * rect2.height) - overlapArea;
    
    return overlapArea / totalArea;
  }

  checkRuleCompliance() {
    const issues = [];
    
    // Check CSS custom properties usage
    const cssVars = this.getCSSVariables();
    if (!cssVars['--safe-area-right']) issues.push('Missing --safe-area-right variable');
    if (!cssVars['--z-sidebar']) issues.push('Missing z-index variables');
    
    // Check safe areas
    const textInput = document.getElementById('object-input');
    if (textInput && !getComputedStyle(textInput).marginRight.includes('var(--safe-area-right)')) {
      issues.push('Text input not using safe area margin');
    }
    
    // Check for borders
    const elements = document.querySelectorAll('button, .mode-label, .btn');
    elements.forEach(el => {
      const borderWidth = getComputedStyle(el).borderWidth;
      if (borderWidth !== '0px') {
        issues.push(`Element ${el.className} has border: ${borderWidth}`);
      }
    });
    
    // Check black backgrounds
    const containers = document.querySelectorAll('body, .app-container, .analysis-section');
    containers.forEach(el => {
      const bgColor = getComputedStyle(el).backgroundColor;
      if (bgColor !== 'rgb(0, 0, 0)' && bgColor !== 'rgba(0, 0, 0, 1)') {
        issues.push(`Element ${el.className} not black background: ${bgColor}`);
      }
    });
    
    return {
      compliant: issues.length === 0,
      issues: issues,
      score: Math.max(0, 100 - (issues.length * 10))
    };
  }

  generatePromptContext() {
    const uiSnapshot = this.generateUISnapshot();
    const ruleCompliance = this.checkRuleCompliance();
    const recentChanges = this.changeHistory.slice(-5);
    
    return {
      current_ui_state: uiSnapshot,
      rule_compliance: ruleCompliance,
      recent_changes: recentChanges,
      active_issues: this.getActiveIssues(uiSnapshot, ruleCompliance),
      suggested_focus: this.suggestNextActions(ruleCompliance, uiSnapshot)
    };
  }

  getActiveIssues(uiSnapshot, ruleCompliance) {
    const issues = [];
    
    // Add rule compliance issues
    ruleCompliance.issues.forEach(issue => {
      issues.push({ type: 'compliance', description: issue, priority: 'high' });
    });
    
    // Add overlap issues
    uiSnapshot.overlappingElements.forEach(overlap => {
      issues.push({ 
        type: 'overlap', 
        description: `${overlap.element1} overlaps ${overlap.element2}`,
        priority: overlap.severity > 0.1 ? 'high' : 'medium'
      });
    });
    
    // Add error state issues
    uiSnapshot.activeElements.errorStates.forEach(error => {
      issues.push({ type: 'error', description: error.message, priority: 'high' });
    });
    
    return issues;
  }

  suggestNextActions(ruleCompliance, uiSnapshot) {
    const suggestions = [];
    
    if (ruleCompliance.score < 90) {
      suggestions.push('Fix rule compliance issues first');
    }
    
    if (uiSnapshot.overlappingElements.length > 0) {
      suggestions.push('Address element overlaps using layout debug mode');
    }
    
    if (uiSnapshot.responsiveState === 'mobile' && uiSnapshot.layout.sidebarVisible) {
      suggestions.push('Test mobile layout with sidebar closed');
    }
    
    return suggestions;
  }

  exportContextForPrompt() {
    const context = this.generatePromptContext();
    const formatted = this.formatContextForPrompt(context);
    
    navigator.clipboard.writeText(formatted).then(() => {
      
    });
    
    return formatted;
  }

  formatContextForPrompt(context) {
    return `
## Current UI State (Auto-Generated ${context.current_ui_state.timestamp})
**Viewport**: ${context.current_ui_state.viewport.width}x${context.current_ui_state.viewport.height} (${context.current_ui_state.responsiveState})
**Rule Compliance**: ${context.rule_compliance.score}%

**Active Issues** (${context.active_issues.length}):
${context.active_issues.map(issue => `- [${issue.priority.toUpperCase()}] ${issue.type}: ${issue.description}`).join('\n')}

**Layout State**:
- Sidebar visible: ${context.current_ui_state.layout.sidebarVisible}
- Safe area right: ${context.current_ui_state.layout.safeAreaRight}
- Focused element: ${context.current_ui_state.activeElements.focusedElement || 'none'}
- Visible sections: ${context.current_ui_state.activeElements.visibleSections.join(', ')}

**Overlapping Elements**: ${context.current_ui_state.overlappingElements.length > 0 ? 'DETECTED' : 'None'}
${context.current_ui_state.overlappingElements.map(overlap => `  - ${overlap.element1} ↔ ${overlap.element2} (severity: ${(overlap.severity * 100).toFixed(1)}%)`).join('\n')}

**Suggested Focus**: ${context.suggested_focus.join(', ')}

**Recent Changes** (${context.recent_changes.length}):
${context.recent_changes.map(change => `- ${change.timestamp}: ${change.description}`).join('\n')}
    `.trim();
  }

  captureUIState(label) {
    const snapshot = {
      label: label,
      timestamp: Date.now(),
      context: this.generatePromptContext()
    };
    
    this.snapshots.set(label, snapshot);
    
    return snapshot;
  }

  setupMutationObserver() {
    new MutationObserver((mutations) => {
      if (this.isSignificantChange(mutations)) {
        this.recordChange('DOM mutation detected');
      }
    }).observe(document.body, { 
      childList: true, 
      subtree: true, 
      attributes: true,
      attributeFilter: ['class', 'style']
    });
  }

  setupEventListeners() {
    window.addEventListener('resize', this.debounce(() => {
      this.recordChange(`Viewport resized to ${window.innerWidth}x${window.innerHeight}`);
    }, 500));
  }

  isSignificantChange(mutations) {
    return mutations.some(mutation => 
      mutation.type === 'childList' && mutation.addedNodes.length > 0 ||
      mutation.type === 'attributes' && ['class', 'style'].includes(mutation.attributeName)
    );
  }

  recordChange(description) {
    this.changeHistory.push({
      timestamp: new Date().toISOString(),
      description: description
    });
    
    // Keep only last 20 changes
    if (this.changeHistory.length > 20) {
      this.changeHistory = this.changeHistory.slice(-20);
    }
  }

  debounce(func, wait) {
    let timeout;
    return function executedFunction(...args) {
      const later = () => {
        clearTimeout(timeout);
        func(...args);
      };
      clearTimeout(timeout);
      timeout = setTimeout(later, wait);
    };
  }
}

// Initialize and expose globally
window.uiContextAutomation = new UIContextAutomation();

// Methods will be added to main app when it initializes 