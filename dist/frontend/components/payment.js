class PaymentManager {
  constructor(mode = 'modal') {
    this.stripe = null;
    this.elements = null;
    this.cardElement = null;
    this.paymentRequest = null;
    this.setupInProgress = false;
    this.mode = mode; // 'modal' or 'sidebar'
  }

  isLocalDevelopment() {
     is deprecated - use developer account instead');
    return false;
  }

  isDeveloperAccount() {
    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    if (!deviceToken || !cardInfo) {
      return false;
    }
    
    try {
      const card = JSON.parse(cardInfo);
      return card.brand === 'Development' || 
             card.last4 === '0000' ||
             deviceToken.includes('local_dev_token') ||
             deviceToken.includes('dev_');
    } catch (error) {
      return false;
    }
  }

  async checkInitialPaymentSetup() {
    // Check for actual payment method regardless of toggle
    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    if (!deviceToken || !cardInfo) {
      
      this.showPaymentModal();
      return false;
    }
    
    
    return true;
  }

  showPaymentModal() {
    if (this.mode === 'sidebar') {
      return this.showPaymentSection();
    }
    
    
    const modal = document.getElementById('payment-modal');
    const backdrop = document.getElementById('modal-backdrop');
    const mainInterface = document.getElementById('main-app-interface');
    
    
    
    if (backdrop) {
      backdrop.classList.remove('hidden');
      setTimeout(() => {
        backdrop.classList.add('show');
      }, 10);
    }
    
    if (modal) {
      modal.style.display = 'block';
      modal.classList.remove('hidden');
      
      setTimeout(() => {
        modal.classList.add('show');
        
      }, 10);
    } else {
      console.error('❌ Payment modal element not found!');
    }
    
    if (mainInterface) {
      mainInterface.classList.add('payment-required');
      mainInterface.classList.add('modal-showing');
    }
    
    this.initializePaymentSetup();
  }

  // Sidebar mode payment section
  showPaymentSection() {
    
    const paymentSection = document.getElementById('payment-section');
    
    if (paymentSection) {
      paymentSection.classList.remove('hidden');
      paymentSection.style.display = 'block';
      
      setTimeout(() => {
        paymentSection.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
      }, 100);
    }
    
    this.initializePaymentSetup();
  }

  hidePaymentSection() {
    
    const paymentSection = document.getElementById('payment-section');
    
    if (paymentSection) {
      paymentSection.classList.add('hidden');
      paymentSection.style.display = 'none';
    }
  }

  // Simplified payment check for sidebar mode
  async checkPaymentRequired() {
    if (this.mode === 'modal') {
      return this.checkInitialPaymentSetup();
    }
    
    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    if (!deviceToken || !cardInfo) {
      
      this.updateAccountButton();
      this.hidePaymentSection();
      return false;
    }
    
    
    this.updateAccountButton();
    return true;
  }

  // Update account button to show payment status (sidebar mode)
  updateAccountButton() {
    const accountStatus = document.getElementById('account-status');
    const accountName = document.getElementById('account-name');
    
    if (accountStatus && accountName) {
      const cardInfo = localStorage.getItem('molCardInfo');
      
      if (cardInfo) {
        try {
          const parsed = JSON.parse(cardInfo);
          accountName.textContent = parsed.name || 'Account Active';
          accountStatus.style.color = '#00d4ff';
          
          if (parsed.usage !== undefined) {
            accountName.textContent = `${parsed.name || 'Active'} (${parsed.usage} analyses)`;
          }
        } catch (e) {
          accountName.textContent = 'Account Active';
          accountStatus.style.color = '#00d4ff';
        }
      } else {
        accountName.textContent = '';
        accountStatus.style.color = '#ff0000';
      }
      
      accountStatus.style.cursor = 'pointer';
      accountStatus.onclick = () => {
        
        const deviceToken = localStorage.getItem('molDeviceToken');
        if (deviceToken) {
          this.showCardManagement();
        } else {
          this.showPaymentSection();
        }
      };
    }
  }

  // Simple card management for sidebar mode
  showCardManagement() {
    if (this.mode === 'modal') {
      return this.showCardManagementModal();
    }
    
    const cardInfo = localStorage.getItem('molCardInfo');
    if (cardInfo) {
      try {
        const parsed = JSON.parse(cardInfo);
        alert(`Account: ${parsed.name || 'Active'}\nUsage: ${parsed.usage || 0} analyses\n\nPayment is active and working.`);
      } catch (e) {
        alert('Payment is active and working.');
      }
    }
  }



  hidePaymentModal() {
    
    const modal = document.getElementById('payment-modal');
    const backdrop = document.getElementById('modal-backdrop');
    const mainInterface = document.getElementById('main-app-interface');
    
    if (backdrop) {
      backdrop.classList.remove('show');
      setTimeout(() => {
        backdrop.classList.add('hidden');
      }, 300);
    }
    
    if (modal) {
      modal.classList.remove('show');
      setTimeout(() => {
        modal.style.display = 'none';
        modal.classList.add('hidden');
      }, 300);
    }
    
    if (mainInterface) {
      mainInterface.classList.remove('payment-required');
      mainInterface.classList.remove('modal-showing');
    }
    
    
  }



  updateAccountStatus(user) {
    const accountStatus = document.getElementById('account-status');
    const accountName = document.getElementById('account-name');
    
    if (accountStatus) {
      // Always show the account status
      accountStatus.classList.remove('hidden');
      accountStatus.style.display = 'flex';
      
      if (user && user.name) {
        accountName.textContent = user.name;
        accountStatus.style.color = '#00d4ff'; // Blue when set up
      } else {
        accountName.textContent = 'Add Card';
        accountStatus.style.color = '#ffa500'; // Orange when not set up
      }
      
      // Add click handler to show payment modal when no payment set up
      const deviceToken = localStorage.getItem('molDeviceToken');
      const cardInfo = localStorage.getItem('molCardInfo');
      
      if (!deviceToken || !cardInfo) {
        // No payment set up - make it clickable to show modal
        accountStatus.style.cursor = 'pointer';
        accountStatus.onclick = () => {
          this.showPaymentModal();
        };
      } else {
        // Payment set up - make it clickable to show card management
        accountStatus.style.cursor = 'pointer';
        accountStatus.onclick = () => {
          this.showCardManagementModal();
        };
      }
    }
  }

  // Card Management Modal Methods
  showCardManagementModal() {
    
    const modal = document.getElementById('card-management-modal');
    const backdrop = document.getElementById('modal-backdrop');
    const mainInterface = document.getElementById('main-app-interface');
    
    if (backdrop) {
      backdrop.classList.remove('hidden');
      setTimeout(() => {
        backdrop.classList.add('show');
      }, 10);
    }
    
    if (modal) {
      modal.style.display = 'block';
      modal.classList.remove('hidden');
      
      setTimeout(() => {
        modal.classList.add('show');
        
      }, 10);
    } else {
      console.error('❌ Card management modal element not found!');
    }
    
    if (mainInterface) {
      mainInterface.classList.add('modal-showing');
    }
    
    this.loadUserCards();
  }



  hideCardManagementModal() {
    
    const modal = document.getElementById('card-management-modal');
    const backdrop = document.getElementById('modal-backdrop');
    const mainInterface = document.getElementById('main-app-interface');
    
    if (backdrop) {
      backdrop.classList.remove('show');
      setTimeout(() => {
        backdrop.classList.add('hidden');
      }, 300);
    }
    
    if (modal) {
      modal.classList.remove('show');
      setTimeout(() => {
        modal.style.display = 'none';
        modal.classList.add('hidden');
      }, 300);
    }
    
    if (mainInterface) {
      mainInterface.classList.remove('modal-showing');
    }
    
    
  }

  async loadUserCards() {
    const deviceToken = localStorage.getItem('molDeviceToken');
    if (!deviceToken) {
      console.error('No device token found');
      return;
    }

    try {
      const response = await fetch(`/get-payment-methods?device_token=${encodeURIComponent(deviceToken)}`, {
        method: 'GET',
        headers: { 'Content-Type': 'application/json' }
      });

      if (!response.ok) {
        throw new Error(`Failed to load cards: ${response.status}`);
      }

      const result = await response.json();
      this.displayCards(result.payment_methods, result.default_method);
      
    } catch (error) {
      console.error('Failed to load user cards:', error);
      this.showCardError('Failed to load payment methods');
    }
  }

  displayCards(cards, defaultMethodId) {
    const cardsList = document.getElementById('cards-list');
    if (!cardsList) return;

    cardsList.innerHTML = '';

    if (!cards || cards.length === 0) {
      cardsList.innerHTML = `
        <div class="empty-state">
          <p>No payment methods</p>
          <button type="button" class="text-button primary" onclick="paymentManager.showAddCardForm()">Add your first payment method</button>
        </div>
      `;
      return;
    }

    cards.forEach(card => {
      const methodElement = document.createElement('div');
      methodElement.className = `payment-method ${card.is_default ? 'default' : ''}`;
      methodElement.innerHTML = `
        <div class="method-info">
          <div class="method-primary">
            <span class="method-brand">${this.formatCardBrand(card.brand)}</span>
            <span class="method-number">•••• ${card.last4}</span>
            ${card.is_default ? '<span class="default-label">Default</span>' : ''}
          </div>
          <div class="method-secondary">
            Expires ${card.exp_month}/${card.exp_year}
          </div>
        </div>
        <div class="method-actions">
          <button type="button" class="action-link" onclick="paymentManager.editCard('${card.id}')">Edit</button>
          ${!card.is_default ? `<button type="button" class="action-link" onclick="paymentManager.setDefaultCard('${card.id}')">Set as default</button>` : ''}
          <button type="button" class="action-link delete" onclick="paymentManager.deleteCard('${card.id}')">Delete</button>
        </div>
      `;
      cardsList.appendChild(methodElement);
    });
  }

  formatCardBrand(brand) {
    const brands = {
      'visa': 'Visa',
      'mastercard': 'Mastercard',
      'amex': 'American Express',
      'discover': 'Discover',
      'diners': 'Diners Club',
      'jcb': 'JCB',
      'unionpay': 'UnionPay'
    };
    return brands[brand] || brand.charAt(0).toUpperCase() + brand.slice(1);
  }

  showAddCardForm() {
    const currentCards = document.getElementById('current-cards');
    const editForm = document.getElementById('edit-card-form');
    const editTitle = document.getElementById('edit-card-title');
    
    if (currentCards) currentCards.classList.add('hidden');
    if (editForm) editForm.classList.remove('hidden');
    if (editTitle) editTitle.textContent = 'Add New Payment Method';
    
    this.currentEditingCardId = null;
    this.setupCardForm('edit-card-element');
    this.setupEditFormSubmission();
  }

  editCard(cardId) {
    const currentCards = document.getElementById('current-cards');
    const editForm = document.getElementById('edit-card-form');
    const editTitle = document.getElementById('edit-card-title');
    
    if (currentCards) currentCards.classList.add('hidden');
    if (editForm) editForm.classList.remove('hidden');
    if (editTitle) editTitle.textContent = 'Edit Payment Method';
    
    this.currentEditingCardId = cardId;
    this.setupCardForm('edit-card-element');
    this.setupEditFormSubmission();
  }

  setupEditFormSubmission() {
    const editForm = document.getElementById('card-edit-form');
    if (editForm) {
      // Remove existing event listeners to prevent duplicates
      editForm.removeEventListener('submit', this.handleEditFormSubmission);
      
      // Add new event listener
      this.handleEditFormSubmission = (event) => {
        event.preventDefault();
        this.saveCardChanges();
      };
      
      editForm.addEventListener('submit', this.handleEditFormSubmission);
      
    }
  }

  async setDefaultCard(cardId) {
    const deviceToken = localStorage.getItem('molDeviceToken');
    if (!deviceToken) return;

    try {
      const response = await fetch('/set-default-payment-method', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ 
          device_token: deviceToken,
          payment_method_id: cardId 
        })
      });

      if (!response.ok) {
        throw new Error(`Failed to set default card: ${response.status}`);
      }

      
      this.loadUserCards(); // Refresh the cards list
      
    } catch (error) {
      console.error('Failed to set default card:', error);
      this.showCardError('Failed to update default card');
    }
  }

  async deleteCard(cardId) {
    if (!confirm('Are you sure you want to delete this payment method?')) {
      return;
    }

    const deviceToken = localStorage.getItem('molDeviceToken');
    if (!deviceToken) return;

    try {
      const response = await fetch('/delete-payment-method', {
        method: 'DELETE',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ 
          device_token: deviceToken,
          payment_method_id: cardId 
        })
      });

      if (!response.ok) {
        throw new Error(`Failed to delete card: ${response.status}`);
      }

      
      this.loadUserCards(); // Refresh the cards list
      
      // If this was the only card, clear local storage and show payment setup
      const result = await response.json();
      if (result.message && result.message.includes('last card')) {
        localStorage.removeItem('molDeviceToken');
        localStorage.removeItem('molCardInfo');
        this.hideCardManagementModal();
        this.showPaymentModal();
      }
      
    } catch (error) {
      console.error('Failed to delete card:', error);
      this.showCardError('Failed to delete payment method');
    }
  }

  cancelCardEdit() {
    const currentCards = document.getElementById('current-cards');
    const editForm = document.getElementById('edit-card-form');
    
    if (currentCards) currentCards.classList.remove('hidden');
    if (editForm) editForm.classList.add('hidden');
    
    // Clean up edit elements
    if (this.editCardElement) {
      this.editCardElement.destroy();
      this.editCardElement = null;
    }
    this.editElements = null;
    
    // Remove form event listener
    const cardEditForm = document.getElementById('card-edit-form');
    if (cardEditForm && this.handleEditFormSubmission) {
      cardEditForm.removeEventListener('submit', this.handleEditFormSubmission);
    }
    
    // Clear any error messages
    const errorElement = document.getElementById('edit-card-errors');
    if (errorElement) {
      errorElement.textContent = '';
      errorElement.style.display = 'none';
    }
    
    this.currentEditingCardId = null;
  }

  setupCardForm(elementId) {
    // Create new Stripe Elements instance for card editing
    if (!this.stripe) {
      console.error('Stripe not initialized for card form');
      return;
    }
    
    const editElements = this.stripe.elements({
      appearance: {
        theme: 'night',
        variables: {
          colorPrimary: '#00d4ff',
          colorBackground: 'rgba(255, 255, 255, 0.08)',
          colorText: '#ffffff',
          colorDanger: '#ff6b6b',
          borderRadius: '8px'
        }
      }
    });
    
    const editCardElement = editElements.create('payment', {
      layout: 'tabs'
    });
    
    const container = document.getElementById(elementId);
    if (container) {
      container.innerHTML = ''; // Clear existing content
      editCardElement.mount(`#${elementId}`);
      
      // Store reference for form submission
      this.editElements = editElements;
      this.editCardElement = editCardElement;
      
      // Set up error handling for edit form
      editCardElement.on('change', (event) => {
        const errorElement = document.getElementById('edit-card-errors');
        if (errorElement) {
          if (event.error) {
            errorElement.textContent = event.error.message;
            errorElement.style.display = 'block';
          } else {
            errorElement.textContent = '';
            errorElement.style.display = 'none';
          }
        }
      });
    }
  }

  async saveCardChanges() {
    if (!this.stripe || !this.editElements) {
      console.error('Stripe not initialized for card editing');
      return;
    }

    const deviceToken = localStorage.getItem('molDeviceToken');
    const userName = document.getElementById('edit-user-name').value;
    
    if (!deviceToken) return;

    try {
      const saveBtn = document.getElementById('save-card-btn');
      const btnText = saveBtn.querySelector('.btn-text');
      const btnLoading = saveBtn.querySelector('.btn-loading');
      
      btnText.classList.add('hidden');
      btnLoading.classList.remove('hidden');
      saveBtn.disabled = true;

      // Trigger form validation and wallet collection
      const { error: submitError } = await this.editElements.submit();
      if (submitError) {
        this.handleEditError(submitError);
        return;
      }

      // Create PaymentMethod
      const { error, paymentMethod } = await this.stripe.createPaymentMethod({
        elements: this.editElements,
        params: {
          billing_details: {
            name: userName || 'Molecular Analysis User',
          }
        }
      });

      if (error) {
        this.handleEditError(error);
        return;
      }

      const endpoint = this.currentEditingCardId ? '/update-payment-method' : '/setup-payment-method';
      const payload = {
        device_token: deviceToken,
        payment_method: paymentMethod.id,
        name: userName
      };

      const response = await fetch(endpoint, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(payload)
      });

      if (!response.ok) {
        throw new Error(`Failed to save card: ${response.status}`);
      }

      const result = await response.json();

      // Handle 3D Secure if required
      if (result.requires_action && result.client_secret) {
        const { error: confirmError } = await this.stripe.confirmCardSetup(
          result.client_secret
        );
        
        if (confirmError) {
          throw new Error(confirmError.message);
        }
      }

      
      this.cancelCardEdit();
      this.loadUserCards();
      
    } catch (error) {
      console.error('Failed to save card:', error);
      this.handleEditError(error);
    } finally {
      const saveBtn = document.getElementById('save-card-btn');
      if (saveBtn) {
        const btnText = saveBtn.querySelector('.btn-text');
        const btnLoading = saveBtn.querySelector('.btn-loading');
        
        btnText.classList.remove('hidden');
        btnLoading.classList.add('hidden');
        saveBtn.disabled = false;
      }
    }
  }

  handleEditError(error) {
    let message = 'An unexpected error occurred.';
    
    if (error.type === 'card_error' || error.type === 'validation_error') {
      message = error.message;
    } else if (error.message) {
      message = error.message;
    }
    
    const errorElement = document.getElementById('edit-card-errors');
    if (errorElement) {
      errorElement.textContent = message;
      errorElement.style.display = 'block';
    }
    
    console.error('Edit card error:', error);
  }

  showCardError(message) {
    // Simple error display - could be enhanced with a toast notification
    alert(message);
  }

  isPaymentRequired() {
    const paymentModal = document.getElementById('payment-modal');
    return paymentModal && paymentModal.classList.contains('show');
  }

  async checkPaymentMethod() {
    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    if (!deviceToken || !cardInfo) {
      this.showPaymentModal();
      return false;
    }
    
    try {
      const response = await fetch('/validate-payment', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ device_token: deviceToken })
      });
      
      if (!response.ok) {
        localStorage.removeItem('molDeviceToken');
        localStorage.removeItem('molCardInfo');
        this.showPaymentModal();
        return false;
      }
      
      const result = await response.json();
      const localCardInfo = JSON.parse(cardInfo);
      localCardInfo.usage = result.user.usage;
      localStorage.setItem('molCardInfo', JSON.stringify(localCardInfo));
      
      return true;
      
    } catch (error) {
      console.error('Payment validation error:', error);
      return true;
    }
  }

  async incrementUsage() {
    const deviceToken = localStorage.getItem('molDeviceToken');
    if (!deviceToken) return;
    
    try {
      const response = await fetch('/increment-usage', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ device_token: deviceToken })
      });
      
      if (response.ok) {
        const result = await response.json();
        const cardInfo = JSON.parse(localStorage.getItem('molCardInfo') || '{}');
        cardInfo.usage = result.usage;
        localStorage.setItem('molCardInfo', JSON.stringify(cardInfo));
        
        
      }
    } catch (error) {
      console.error('Usage increment error:', error);
      const cardInfo = JSON.parse(localStorage.getItem('molCardInfo') || '{}');
      cardInfo.usage = (cardInfo.usage || 0) + 1;
      localStorage.setItem('molCardInfo', JSON.stringify(cardInfo));
    }
  }

  resetLocalDevUser() {
     is deprecated - use developer account instead');
    this.setupDeveloperAccount();
  }

  async initializePaymentSetup() {
    
    
    try {
      // Get Stripe publishable key from server
      const response = await fetch('/stripe-config');
      if (!response.ok) {
        throw new Error('Failed to get Stripe configuration');
      }
      
      const config = await response.json();
      
      // Initialize Stripe
      this.stripe = Stripe(config.publishableKey);
      
      
      // Create Elements instance
      this.elements = this.stripe.elements({
        appearance: {
          theme: 'night',
          variables: {
            colorPrimary: '#00d4ff',
            colorBackground: 'rgba(255, 255, 255, 0.08)',
            colorText: '#ffffff',
            colorDanger: '#ff6b6b',
            borderRadius: '8px'
          }
        }
      });
      
      // Create and mount Card Element
      this.cardElement = this.elements.create('card', {
        style: {
          base: {
            color: '#ffffff',
            fontSize: '16px',
            '::placeholder': {
              color: 'rgba(255, 255, 255, 0.7)'
            }
          }
        }
      });
      
      const cardElementContainer = document.getElementById('card-element');
      if (cardElementContainer) {
        this.cardElement.mount('#card-element');
        
      }
      
      // Set up error handling
      this.cardElement.on('change', (event) => {
        this.handleElementChange(event);
      });
      
      // Set up form submission
      this.setupFormSubmission();
      
    } catch (error) {
      console.error('❌ Failed to initialize Stripe:', error);
      this.showError('Payment system unavailable. Please try again.');
    }
  }

  handleElementChange(event) {
    const errorElement = document.getElementById('card-errors');
    if (event.error) {
      errorElement.textContent = event.error.message;
      errorElement.style.display = 'block';
    } else {
      errorElement.textContent = '';
      errorElement.style.display = 'none';
    }
  }

  setupFormSubmission() {
    const form = document.getElementById('card-setup-form');
    if (form) {
      form.addEventListener('submit', (event) => {
        this.handleFormSubmission(event);
      });
      
    }
  }

  async handleFormSubmission(event) {
    event.preventDefault();
    
    if (!this.stripe || !this.elements) {
      console.error('Stripe not initialized');
      return;
    }
    
    const setupBtn = document.getElementById('setup-btn');
    const btnText = setupBtn.querySelector('.btn-text');
    const btnLoading = setupBtn.querySelector('.btn-loading');
    
    // Show loading state
    btnText.classList.add('hidden');
    btnLoading.classList.remove('hidden');
    setupBtn.disabled = true;
    
    try {
      // Trigger form validation and wallet collection
      const { error: submitError } = await this.elements.submit();
      if (submitError) {
        this.handleError(submitError);
        return;
      }
      
      // Create PaymentMethod
      const { error, paymentMethod } = await this.stripe.createPaymentMethod({
        elements: this.elements,
        params: {
          billing_details: {
            name: document.getElementById('user-name').value || 'Molecular Analysis User',
          }
        }
      });
      
      if (error) {
        this.handleError(error);
        return;
      }
      
      
      
      // Send to server for setup
      await this.setupPaymentMethodOnServer(paymentMethod);
      
    } catch (error) {
      console.error('Payment setup error:', error);
      this.handleError(error);
    } finally {
      // Reset button state
      btnText.classList.remove('hidden');
      btnLoading.classList.add('hidden');
      setupBtn.disabled = false;
    }
  }

  async setupPaymentMethodOnServer(paymentMethod) {
    try {
      const response = await fetch('/setup-payment-method', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({
          payment_method: paymentMethod.id,
          name: document.getElementById('user-name').value || null
        })
      });
      
      if (!response.ok) {
        throw new Error(`Server error: ${response.status}`);
      }
      
      const result = await response.json();
      
      // Handle 3D Secure if required
      if (result.requires_action && result.client_secret) {
        const { error: confirmError } = await this.stripe.confirmCardSetup(
          result.client_secret
        );
        
        if (confirmError) {
          throw new Error(confirmError.message);
        }
      }
      
      // Store payment info locally
      this.storePaymentInfo(paymentMethod);
      
      // Show success
      this.showSuccess();
      
    } catch (error) {
      throw new Error(`Failed to setup payment method: ${error.message}`);
    }
  }

  storePaymentInfo(paymentMethod) {
    const deviceToken = this.generateDeviceToken();
    const cardInfo = {
      last4: paymentMethod.card ? paymentMethod.card.last4 : '••••',
      brand: paymentMethod.card ? paymentMethod.card.brand : 'card',
      usage: 0,
      name: document.getElementById('user-name').value || null,
      setupDate: new Date().toISOString()
    };
    
    localStorage.setItem('molDeviceToken', deviceToken);
    localStorage.setItem('molCardInfo', JSON.stringify(cardInfo));
    
    
  }

  generateDeviceToken() {
    return 'mol_' + Date.now() + '_' + Math.random().toString(36).substr(2, 9);
  }

  showSuccess() {
    const form = document.getElementById('card-setup-form');
    const success = document.getElementById('setup-success');
    
    if (form && success) {
      form.classList.add('hidden');
      success.classList.remove('hidden');
      
      // Update account status
      const cardInfo = JSON.parse(localStorage.getItem('molCardInfo'));
      this.updateAccountStatus({ name: cardInfo.name });
      
      
    }
  }

  handleError(error) {
    let message = 'An unexpected error occurred.';
    
    if (error.type === 'card_error' || error.type === 'validation_error') {
      message = error.message;
    } else if (error.message) {
      message = error.message;
    }
    
    const errorElement = document.getElementById('card-errors');
    if (errorElement) {
      errorElement.textContent = message;
      errorElement.style.display = 'block';
    }
    
    console.error('Payment error:', error);
  }

  showError(message) {
    const errorElement = document.getElementById('card-errors');
    if (errorElement) {
      errorElement.textContent = message;
      errorElement.style.display = 'block';
    }
  }

  setupLocalDevUser() {
     is deprecated - use developer account instead');
    this.setupDeveloperAccount();
  }

  getLocalDevUser() {
     is deprecated - use developer account instead');
    return this.getDeveloperAccount();
  }

  setupDeveloperAccount() {
    const developerUser = {
      name: 'Developer',
      email: 'developer@mol.com',
      usage: 0,
      device_token: 'dev_token_' + Date.now(),
      card_info: {
        last4: '0000',
        brand: 'Development',
        usage: 0
      }
    };

    localStorage.setItem('molDeviceToken', developerUser.device_token);
    localStorage.setItem('molCardInfo', JSON.stringify(developerUser.card_info));
    localStorage.setItem('molDeveloperUser', JSON.stringify(developerUser));

    this.updateAccountStatus(developerUser);
    
    
    return developerUser;
  }

  getDeveloperAccount() {
    const developerUser = localStorage.getItem('molDeveloperUser');
    return developerUser ? JSON.parse(developerUser) : null;
  }
  
  // Function to clear payment setup for testing
  clearPaymentSetup() {
    localStorage.removeItem('molDeviceToken');
    localStorage.removeItem('molCardInfo');
    localStorage.removeItem('molDeveloperUser');
    
  }
}

// Export singleton instance
export const paymentManager = new PaymentManager(); 