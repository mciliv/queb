// Fix for multiple enter key triggers in app.js
// Replace the setupTextAnalysis method with this version:

setupTextAnalysis() {
  const textInput = document.getElementById('object-input');
  if (textInput) {
    // Prevent multiple rapid triggers with debouncing
    let lastTriggerTime = 0;
    const DEBOUNCE_MS = 1000; // 1 second cooldown
    
    textInput.addEventListener('keydown', (event) => {
      if (event.key === 'Enter') {
        event.preventDefault();
        
        const now = Date.now();
        if (now - lastTriggerTime < DEBOUNCE_MS) {
          console.log('🚫 Enter key debounced - too soon since last trigger');
          return;
        }
        
        lastTriggerTime = now;
        this.handleTextAnalysis();
      }
    });
  }
}
