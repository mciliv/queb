import React, { useState, useEffect, useRef } from 'react';

const AccountSetup = () => {
    const [stripe, setStripe] = useState(null);
    const [cardElement, setCardElement] = useState(null);
    const [currentUser, setCurrentUser] = useState(null);
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState('');
    const [success, setSuccess] = useState('');
    const cardElementRef = useRef(null);

    useEffect(() => {
        initializeStripe();
        checkExistingUser();
    }, []);

    const initializeStripe = async () => {
        try {
            const response = await fetch('/api/stripe-config');
            const config = await response.json();
            
            const stripeInstance = window.Stripe(config.publishableKey);
            setStripe(stripeInstance);
            
            const elements = stripeInstance.elements();
            const cardEl = elements.create('card', {
                style: {
                    base: {
                        fontSize: '16px',
                        color: '#ffffff',
                        fontFamily: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif',
                        '::placeholder': {
                            color: 'rgba(255, 255, 255, 0.6)'
                        }
                    },
                    invalid: {
                        color: '#ff6b6b',
                        iconColor: '#ff6b6b'
                    }
                },
                hidePostalCode: true
            });
            
            cardEl.mount(cardElementRef.current);
            setCardElement(cardEl);
            
            cardEl.on('change', (event) => {
                setError(event.error ? event.error.message : '');
            });
        } catch (err) {
            setError('Failed to initialize payment system');
        }
    };

    const checkExistingUser = async () => {
        try {
            const response = await fetch('/api/user/check');
            const user = await response.json();
            setCurrentUser(user);
        } catch (err) {
            console.log('No existing user found');
        }
    };

    const handleSubmit = async (event) => {
        event.preventDefault();
        setLoading(true);
        setError('');

        if (!stripe || !cardElement) {
            setError('Payment system not loaded');
            setLoading(false);
            return;
        }

        const { error: stripeError, paymentMethod } = await stripe.createPaymentMethod({
            type: 'card',
            card: cardElement,
        });

        if (stripeError) {
            setError(stripeError.message);
            setLoading(false);
            return;
        }

        try {
            const response = await fetch('/api/create-subscription', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({
                    payment_method: paymentMethod.id,
                }),
            });

            const result = await response.json();

            if (result.error) {
                setError(result.error);
            } else {
                setSuccess('Account created successfully!');
                setCurrentUser(result.user);
            }
        } catch (err) {
            setError('Failed to create account');
        }

        setLoading(false);
    };

    if (currentUser) {
        return (
            <div className="account-container">
                <div className="account-card">
                    <h2>Welcome back!</h2>
                    <p>Your account is active and ready to use.</p>
                    <div className="user-info">
                        <p>Email: {currentUser.email}</p>
                        <p>Status: {currentUser.status}</p>
                    </div>
                </div>
            </div>
        );
    }

    return (
        <div className="account-container">
            <div className="account-card">
                <h2>Set Up Your Account</h2>
                <p>Get unlimited molecular analysis with your subscription</p>
                
                <form onSubmit={handleSubmit}>
                    <div className="payment-section">
                        <label htmlFor="card-element">Card Information</label>
                        <div ref={cardElementRef} id="card-element"></div>
                    </div>
                    
                    {error && <div className="error-message">{error}</div>}
                    {success && <div className="success-message">{success}</div>}
                    
                    <button 
                        type="submit" 
                        disabled={loading || !stripe}
                        className="submit-button"
                    >
                        {loading ? 'Processing...' : 'Start Subscription'}
                    </button>
                </form>
                
                <div className="security-note">
                    <p>🔒 Payments secured by Stripe</p>
                    <p>Cancel anytime from your account</p>
                </div>
            </div>
        </div>
    );
};

export default AccountSetup;