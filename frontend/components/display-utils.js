// Display State Management Utility
// Replaces inline style.display assignments with CSS classes

export class DisplayUtils {
  
  // Show element
  static show(element) {
    if (!element) return;
    element.classList.remove('display-none', 'modal-closed');
    element.classList.add('display-block', 'modal-open');
  }
  
  // Hide element
  static hide(element) {
    if (!element) return;
    element.classList.remove('display-block', 'modal-open');
    element.classList.add('display-none', 'modal-closed');
  }
  
  // Show as flex
  static showFlex(element) {
    if (!element) return;
    element.classList.remove('display-none', 'account-status-hidden');
    element.classList.add('display-flex', 'account-status-visible');
  }
  
  // Toggle visibility
  static toggle(element, show) {
    if (show) {
      this.show(element);
    } else {
      this.hide(element);
    }
  }
  
  // Set payment status color
  static setPaymentStatus(element, status) {
    if (!element) return;
    element.classList.remove('payment-status-active', 'payment-status-inactive', 'payment-status-error');
    element.classList.add(`payment-status-${status}`);
  }
  
  // Set cursor pointer
  static setCursorPointer(element) {
    if (!element) return;
    element.classList.add('cursor-pointer');
  }
  
  // Set loading state for buttons
  static setButtonLoading(button, loading) {
    if (!button) return;
    if (loading) {
      button.classList.add('btn-loading-state');
      button.classList.remove('btn-normal-state');
    } else {
      button.classList.remove('btn-loading-state');
      button.classList.add('btn-normal-state');
    }
  }
  
  // Show/hide error element
  static showError(element, show = true) {
    if (!element) return;
    if (show) {
      element.classList.remove('display-none');
      element.classList.add('display-block');
    } else {
      element.classList.remove('display-block');
      element.classList.add('display-none');
    }
  }
  
  // Modal management
  static openModal(modal) {
    if (!modal) return;
    modal.classList.remove('modal-closed');
    modal.classList.add('modal-open');
  }
  
  static closeModal(modal) {
    if (!modal) return;
    modal.classList.remove('modal-open');
    modal.classList.add('modal-closed');
  }
}

// Global utility functions for backward compatibility
window.showElement = (element) => DisplayUtils.show(element);
window.hideElement = (element) => DisplayUtils.hide(element);
window.toggleElement = (element, show) => DisplayUtils.toggle(element, show); 