require('../../../src/client/ui-context-automation.js');

const UIContextAutomation = window.uiContextAutomation.constructor;

describe('UIContextAutomation core logic', () => {
  let ui;

  beforeEach(() => {
    ui = new UIContextAutomation();
  });

  test('getResponsiveBreakpoint returns correct breakpoint', () => {
    window.innerWidth = 400;
    expect(ui.getResponsiveBreakpoint()).toBe('mobile');
    window.innerWidth = 600;
    expect(ui.getResponsiveBreakpoint()).toBe('tablet');
    window.innerWidth = 800;
    expect(ui.getResponsiveBreakpoint()).toBe('desktop');
  });

  test('isSignificantChange detects added nodes and attribute changes', () => {
    const addMutation = { type: 'childList', addedNodes: [1], attributeName: null };
    expect(ui.isSignificantChange([addMutation])).toBe(true);
    const attrMutation = { type: 'attributes', addedNodes: [], attributeName: 'class' };
    expect(ui.isSignificantChange([attrMutation])).toBe(true);
    const other = { type: 'attributes', addedNodes: [], attributeName: 'id' };
    expect(ui.isSignificantChange([other])).toBe(false);
  });

  test('recordChange stores descriptions and limits history to 20 entries', () => {
    for (let i = 0; i < 25; i++) {
      ui.recordChange(`change ${i}`);
    }
    expect(ui.changeHistory).toHaveLength(20);
    expect(ui.changeHistory[0].description).toBe('change 5');
  });

  describe('debounce', () => {
    beforeAll(() => jest.useFakeTimers());
    afterAll(() => jest.useRealTimers());

    test('debounce delays execution', () => {
      const fn = jest.fn();
      const debounced = ui.debounce(fn, 100);
      debounced();
      expect(fn).not.toHaveBeenCalled();
      jest.advanceTimersByTime(50);
      expect(fn).not.toHaveBeenCalled();
      jest.advanceTimersByTime(50);
      expect(fn).toHaveBeenCalledTimes(1);
    });
  });

  test('formatContextForPrompt produces expected markdown format', () => {
    const context = {
      current_ui_state: {
        timestamp: '2022-01-01T00:00:00Z',
        viewport: { width: 100, height: 200 },
        responsiveState: 'desktop',
        layout: { sidebarVisible: true, safeAreaRight: '20px' },
        activeElements: { focusedElement: 'input', visibleSections: ['sec1'] },
        overlappingElements: [{ element1: 'a', element2: 'b', severity: 0.1 }]
      },
      rule_compliance: { score: 80 },
      active_issues: [{ priority: 'high', type: 'test', description: 'desc' }],
      suggested_focus: ['focus1'],
      recent_changes: [{ timestamp: 'ts', description: 'd' }]
    };
    const markdown = ui.formatContextForPrompt(context);
    expect(markdown).toContain('## Current UI State (Auto-Generated 2022-01-01T00:00:00Z)');
    expect(markdown).toContain('**Viewport**: 100x200 (desktop)');
    expect(markdown).toContain('- [HIGH] test: desc');
    expect(markdown).toContain('Suggested Focus');
    expect(markdown).toContain('Recent Changes');
  });
});
