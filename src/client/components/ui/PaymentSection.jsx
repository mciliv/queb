import React, { useState, useEffect } from 'react';
import { usePayment } from './PaymentContext';

const PaymentSection = () => {
  const { 
    config, 
    isPaymentSetup, 
    setIsPaymentSetup,
    isLoading,
    setIsLoading 
  } = usePayment();
  
  const [cardElement, setCardElement] = useState(null);
  const [userName, setUserName] = useState('');
  const [error, setError] = useState('');

  useEffect(() => {
    if (config.enabled && window.Stripe) {
      // Initialize Stripe
      const stripe = window.Stripe(config.publishableKey);
      const elements = stripe.elements();
      const card = elements.create('card', {
        style: {
          base: {
            fontSize: '16px',
            color: '#1a1a1a',
            '::placeholder': {
              color: '#999'
            }
          }
        }
      });
      card.mount('#card-element');
      setCardElement({ stripe, card });
    }
  }, [config.enabled, config.publishableKey]);

  const handleSubmit = async (e) => {
    e.preventDefault();
    if (!config.enabled || isLoading || !cardElement) return;

    setIsLoading(true);
    setError('');

    try {
      // Create payment method with Stripe
      const { stripe, card } = cardElement;
      const { paymentMethod, error: stripeError } = await stripe.createPaymentMethod({
        type: 'card',
        card: card,
      });

      if (stripeError) {
        setError(stripeError.message);
        return;
      }

      // Get device info for token generation
      const deviceInfo = navigator.userAgent;

      // Setup payment method with backend
      const response = await fetch('/api/setup-payment-method', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          payment_method: paymentMethod.id,
          device_info: deviceInfo,
          name: userName
        })
      });

      const data = await response.json();

      if (!response.ok) {
        setError(data.error || 'Payment setup failed');
        return;
      }

      // Store device token for future requests
      if (data.device_token) {
        localStorage.setItem('device_token', data.device_token);
      }

      setIsPaymentSetup(true);
    } catch (err) {
      setError('Payment setup failed. Please try again.');
    } finally {
      setIsLoading(false);
    }
  };

  if (!config.enabled) {
    return null;
  }

  if (isPaymentSetup) {
    return (
      <div className="payment-section">
        <div className="payment-content">
          <div className="account-status">
            <div className="account-link">
              <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                <rect x="2" y="4" width="20" height="16" rx="2"/>
                <line x1="2" y1="10" x2="22" y2="10"/>
              </svg>
              <span>{userName || 'Account'}</span>
            </div>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div className="payment-section">
      <div className="payment-content">
        <div className="payment-header">
          <h2>Payment Required</h2>
          <p>$0.25 per analysis</p>
        </div>
        
        <form className="payment-form" onSubmit={handleSubmit}>
          <div className="form-row">
            <label>Card</label>
            <div id="card-element" className="card-input"></div>
            {error && <div className="error-text">{error}</div>}
          </div>
          
          <div className="form-row">
            <label>Name</label>
            <input 
              type="text" 
              placeholder="Optional"
              value={userName}
              onChange={(e) => setUserName(e.target.value)}
              disabled={isLoading}
            />
          </div>
          
          <button type="submit" className="submit-btn" disabled={isLoading}>
            <span className={isLoading ? 'hidden' : ''}>Setup</span>
            <span className={isLoading ? '' : 'hidden'}>...</span>
          </button>
        </form>
      </div>
    </div>
  );
};

export default PaymentSection;