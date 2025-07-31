class SimplePaymentManager {
  constructor() {
    this.stripe = null;
    this.elements = null;
    this.cardElement = null;
  }

  // Show payment section at top
  showPaymentSection() {
    console.log('💳 Showing payment section');
    const paymentSection = document.getElementById('payment-section');
    const separatorLine = document.getElementById('separator-line');
    
    if (paymentSection) {
      paymentSection.classList.remove('hidden');
    }
    
    if (separatorLine) {
      separatorLine.style.display = 'block';
    }
    
    this.initializeStripe();
  }

  // Hide payment section
  hidePaymentSection() {
    console.log('💳 Hiding payment section');
    const paymentSection = document.getElementById('payment-section');
    const separatorLine = document.getElementById('separator-line');
    
    if (paymentSection) {
      paymentSection.classList.add('hidden');
    }
    
    if (separatorLine) {
      separatorLine.style.display = 'none';
    }
  }

  // Check if payment is needed
  async checkPaymentRequired() {
    // Check if payment is globally disabled
    const paymentEnabled = window.app && window.app.paymentEnabled;
    if (!paymentEnabled) {
      console.log('💳 Payment disabled globally - hiding payment section');
      this.hidePaymentSection();
      return true;
    }
    
    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    if (!deviceToken || !cardInfo) {
      console.log('🔧 No payment method found - showing payment section');
      this.showPaymentSection();
      return false;
    }
    
    console.log('✅ Payment method found');
    this.updateAccountButton();
    return true;
  }

  // Initialize Stripe
  async initializeStripe() {
    if (this.stripe) return; // Already initialized
    
    try {
      // Get Stripe config from server
      const response = await fetch('/stripe-config');
      if (!response.ok) {
        throw new Error('Failed to get Stripe config');
      }
      
      const config = await response.json();
      this.stripe = Stripe(config.publishableKey);
      
      // Create elements
      this.elements = this.stripe.elements({
        appearance: {
          theme: 'night',
          variables: {
            colorPrimary: '#00d4ff',
            colorBackground: 'rgba(255, 255, 255, 0.1)',
            colorText: '#ffffff',
            borderRadius: '8px'
          }
        }
      });
      
      // Create and mount card element
      this.cardElement = this.elements.create('payment', {
        layout: 'tabs'
      });
      
      this.cardElement.mount('#card-element');
      
      // Handle errors
      this.cardElement.on('change', (event) => {
        const errorElement = document.getElementById('card-errors');
        if (event.error) {
          errorElement.textContent = event.error.message;
        } else {
          errorElement.textContent = '';
        }
      });
      
      // Handle form submission
      this.setupFormHandler();
      
      console.log('✅ Stripe initialized');
      
    } catch (error) {
      console.error('❌ Stripe initialization failed:', error);
      const errorElement = document.getElementById('card-errors');
      if (errorElement) {
        errorElement.textContent = 'Payment system unavailable. Please try again.';
      }
    }
  }

  // Setup form submission
  setupFormHandler() {
    const form = document.getElementById('payment-form');
    if (!form) return;
    
    form.addEventListener('submit', async (event) => {
      event.preventDefault();
      await this.handlePaymentSetup();
    });
    
    // Close button
    const closeBtn = document.getElementById('payment-close');
    if (closeBtn) {
      closeBtn.addEventListener('click', () => {
        this.hidePaymentSection();
      });
    }
  }

  // Handle payment setup
  async handlePaymentSetup() {
    if (!this.stripe || !this.elements) {
      console.error('Stripe not initialized');
      return;
    }
    
    const submitBtn = document.getElementById('setup-btn');
    const btnText = submitBtn.querySelector('.btn-text');
    const btnLoading = submitBtn.querySelector('.btn-loading');
    const errorElement = document.getElementById('card-errors');
    
    // Show loading
    btnText.classList.add('hidden');
    btnLoading.classList.remove('hidden');
    submitBtn.disabled = true;
    errorElement.textContent = '';
    
    try {
      // Submit elements
      const { error: submitError } = await this.elements.submit();
      if (submitError) {
        throw submitError;
      }
      
      // Create payment method
      const { error, paymentMethod } = await this.stripe.createPaymentMethod({
        elements: this.elements,
        params: {
          billing_details: {
            name: document.getElementById('user-name').value || 'Molecular Analysis User',
          }
        }
      });
      
      if (error) {
        throw error;
      }
      
      // Send to server
      await this.setupPaymentOnServer(paymentMethod);
      
      // Success - hide payment section
      this.hidePaymentSection();
      this.updateAccountButton();
      
      console.log('✅ Payment setup complete');
      
    } catch (error) {
      console.error('❌ Payment setup failed:', error);
      this.showError(error.message || 'Payment setup failed. Please try again.');
    } finally {
      // Reset button
      btnText.classList.remove('hidden');
      btnLoading.classList.add('hidden');
      submitBtn.disabled = false;
    }
  }

  // Send payment method to server
  async setupPaymentOnServer(paymentMethod) {
    const response = await fetch('/setup-payment-method', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        payment_method: paymentMethod.id,
        name: document.getElementById('user-name').value || null
      })
    });
    
    if (!response.ok) {
      throw new Error(`Server error: ${response.status}`);
    }
    
    const result = await response.json();
    
    // Handle 3D Secure if needed
    if (result.requires_action && result.client_secret) {
      const { error: confirmError } = await this.stripe.confirmCardSetup(
        result.client_secret
      );
      
      if (confirmError) {
        throw confirmError;
      }
    }
    
    // Store payment info
    this.storePaymentInfo(paymentMethod);
  }

  // Store payment info locally
  storePaymentInfo(paymentMethod) {
    const deviceToken = 'mol_' + Date.now() + '_' + Math.random().toString(36).substr(2, 9);
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

  // Update account button
  updateAccountButton() {
    const accountBtn = document.getElementById('account-btn');
    const accountText = document.getElementById('account-text');
    
    if (!accountBtn || !accountText) return;
    
    const cardInfo = localStorage.getItem('molCardInfo');
    
    if (cardInfo) {
      const card = JSON.parse(cardInfo);
      accountText.textContent = card.name || 'Payment Set';
      accountBtn.style.color = '#00ff88'; // Green when setup
    } else {
      accountText.textContent = 'Add Card';
      accountBtn.style.color = '#ffa500'; // Orange when not setup
      
      // Click to show payment
      accountBtn.onclick = () => {
        this.showPaymentSection();
      };
    }
  }

  // Show error message
  showError(message) {
    const errorElement = document.getElementById('card-errors');
    if (errorElement) {
      errorElement.textContent = message;
    }
  }

  // Setup developer account (for testing)
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
    
    this.updateAccountButton();
    this.hidePaymentSection();
    
    console.log('🔧 Developer account setup complete');
    return developerUser;
  }

  // Check if payment method is valid
  async checkPaymentMethod() {
    const deviceToken = localStorage.getItem('molDeviceToken');
    const cardInfo = localStorage.getItem('molCardInfo');
    
    if (!deviceToken || !cardInfo) {
      this.showPaymentSection();
      return false;
    }
    
    // Could add server validation here
    return true;
  }

  // Increment usage counter
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
        
        console.log(`📊 Analysis complete - Total usage: ${result.usage}`);
      }
    } catch (error) {
      console.error('Usage increment error:', error);
      // Update locally as fallback
      const cardInfo = JSON.parse(localStorage.getItem('molCardInfo') || '{}');
      cardInfo.usage = (cardInfo.usage || 0) + 1;
      localStorage.setItem('molCardInfo', JSON.stringify(cardInfo));
    }
  }
}

// Export singleton
export const simplePaymentManager = new SimplePaymentManager(); 