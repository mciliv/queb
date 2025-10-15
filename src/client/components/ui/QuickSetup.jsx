import React, { useState, useEffect } from 'react';
import { useStripeSetup } from '../hooks/useStripeSetup';

const QuickSetup = () => {
    const {
        stripe,
        cardElement,
        currentUser,
        loading,
        error,
        success,
        cardElementRef,
        setLoading,
        setError,
        setSuccess,
        setCurrentUser,
    } = useStripeSetup();

    const [userName, setUserName] = useState('');
    const [showExistingUser, setShowExistingUser] = useState(!!currentUser);

    useEffect(() => {
        setShowExistingUser(!!currentUser);
    }, [currentUser]);

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
            billing_details: {
                name: userName || undefined,
            },
        });

        if (stripeError) {
            setError(stripeError.message);
            setLoading(false);
            return;
        }

        try {
            const response = await fetch('/api/setup-payment', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({
                    payment_method: paymentMethod.id,
                    name: userName,
                }),
            });

            const result = await response.json();

            if (result.error) {
                setError(result.error);
            } else {
                setSuccess(true);
            }
        } catch (err) {
            setError('Failed to setup payment');
        }

        setLoading(false);
    };

    const handleContinue = () => {
        // Navigate back to analysis
        window.location.href = '/';
    };

    const handleUseDifferentCard = () => {
        setShowExistingUser(false);
        setCurrentUser(null);
    };

    if (success) {
        return (
            <div className="container">
                <header>
                    <h1>Quick Setup</h1>
                    <p className="subtitle">Just add your card to start analyzing molecules</p>
                    <a href="/" className="back-link">← Back to Analysis</a>
                </header>

                <div className="setup-container">
                    <div className="setup-success">
                        <div className="success-animation">
                            <svg width="64" height="64" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                                <path d="M9 12l2 2 4-4"/>
                                <circle cx="12" cy="12" r="10"/>
                            </svg>
                        </div>
                        <h3>You're all set!</h3>
                        <p>Start analyzing molecules right away</p>
                        <button onClick={() => window.location.href = '/'} className="btn-primary">Start Analyzing</button>
                    </div>
                </div>
            </div>
        );
    }

    return (
        <div className="container">
            <header>
                <h1>Quick Setup</h1>
                <p className="subtitle">Just add your card to start analyzing molecules</p>
                <a href="/" className="back-link">← Back to Analysis</a>
            </header>

            <div className="setup-container">
                {/* Existing User Recognition */}
                {showExistingUser && currentUser && (
                    <div className="existing-user">
                        <div className="user-card">
                            <svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                                <path d="M20 21v-2a4 4 0 0 0-4-4H8a4 4 0 0 0-4 4v2"/>
                                <circle cx="12" cy="7" r="4"/>
                            </svg>
                            <div className="user-info">
                                <h3>Welcome back!</h3>
                                <p>Card ending in <span>{currentUser.cardLast4 || '••••'}</span></p>
                                <p className="usage">Analyses used: <span>{currentUser.usage || 0}</span></p>
                            </div>
                            <button className="btn-primary" onClick={handleContinue}>Continue</button>
                        </div>
                        <button className="btn-link" onClick={handleUseDifferentCard}>Use different card</button>
                    </div>
                )}

                {/* New User Setup */}
                {!showExistingUser && (
                    <div className="new-user">
                        <div className="setup-info">
                            <div className="price-display">
                                <span className="price">$0.25</span>
                                <span className="per">per analysis</span>
                            </div>
                            <p className="security-note">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                                    <rect x="3" y="11" width="18" height="11" rx="2" ry="2"/>
                                    <circle cx="12" cy="16" r="1"/>
                                    <path d="M7 11V7a5 5 0 0 1 10 0v4"/>
                                </svg>
                                Secured with 3D Secure & encrypted storage
                            </p>
                        </div>

                        {/* Express Payment Options */}
                        <div className="express-payment">
                            <div className="express-title">Quick Payment</div>
                            
                            {/* Apple Pay / Google Pay */}
                            <div id="payment-request-button" className="payment-request-btn"></div>
                            
                            {/* PayPal Express */}
                            <div id="paypal-button-container" className="paypal-container"></div>
                            
                            <div className="or-divider">
                                <span>or</span>
                            </div>
                        </div>
                    
                        <form onSubmit={handleSubmit} className="setup-form">
                            {/* Credit Card with Browser Autofill Support */}
                            <div className="form-group">
                                <label htmlFor="card-element">Card Details</label>
                                
                                {/* Hidden fields for browser autofill detection */}
                                <div className="autofill-detector">
                                    <input type="text" name="cc-name" autoComplete="cc-name" tabIndex="-1" style={{ position: 'absolute', left: '-9999px' }} />
                                    <input type="text" name="cc-number" autoComplete="cc-number" tabIndex="-1" style={{ position: 'absolute', left: '-9999px' }} />
                                    <input type="text" name="cc-exp-month" autoComplete="cc-exp-month" tabIndex="-1" style={{ position: 'absolute', left: '-9999px' }} />
                                    <input type="text" name="cc-exp-year" autoComplete="cc-exp-year" tabIndex="-1" style={{ position: 'absolute', left: '-9999px' }} />
                                    <input type="text" name="cc-csc" autoComplete="cc-csc" tabIndex="-1" style={{ position: 'absolute', left: '-9999px' }} />
                                </div>
                                
                                <div ref={cardElementRef} className="card-element"></div>
                                {error && <div className="error-message" role="alert">{error}</div>}
                                
                                {/* Autofill suggestions */}
                                <div className="autofill-note">
                                    <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                                        <path d="M9 12l2 2 4-4"/>
                                        <circle cx="12" cy="12" r="10"/>
                                    </svg>
                                    <span>Browser autofill supported</span>
                                </div>
                            </div>

                            {/* Optional name for personalization */}
                            <div className="form-group optional">
                                <label htmlFor="user-name">Name (optional)</label>
                                <input 
                                    type="text" 
                                    id="user-name" 
                                    name="name" 
                                    placeholder="For personalized experience" 
                                    autoComplete="name"
                                    value={userName}
                                    onChange={(e) => setUserName(e.target.value)}
                                />
                            </div>

                            <button type="submit" className="btn-primary setup-btn" disabled={loading}>
                                <span className="btn-icon">
                                    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                                        <path d="M12 2l3.09 6.26L22 9.27l-5 4.87 1.18 6.88L12 17.77l-6.18 3.25L7 14.14 2 9.27l6.91-1.01L12 2z"/>
                                    </svg>
                                </span>
                                <span className={loading ? 'btn-text-hidden' : 'btn-text'}>Start Analyzing</span>
                                <span className={loading ? 'btn-loading' : 'btn-loading-hidden'}>
                                    <svg className="spinner" width="16" height="16" viewBox="0 0 24 24">
                                        <circle cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="2" fill="none" strokeDasharray="31.416" strokeDashoffset="31.416">
                                            <animate attributeName="stroke-dasharray" dur="2s" values="0 31.416;15.708 15.708;0 31.416" repeatCount="indefinite"/>
                                            <animate attributeName="stroke-dashoffset" dur="2s" values="0;-15.708;-31.416" repeatCount="indefinite"/>
                                        </circle>
                                    </svg>
                                    Securing...
                                </span>
                            </button>
                        </form>

                        <div className="features">
                            <div className="feature">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                                    <path d="M9 12l2 2 4-4"/>
                                    <circle cx="12" cy="12" r="10"/>
                                </svg>
                                <span>Works with all browser autofills</span>
                            </div>
                            <div className="feature">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                                    <circle cx="12" cy="12" r="3"/>
                                    <path d="M12 1v6M12 17v6M4.22 4.22l4.24 4.24M15.54 15.54l4.24 4.24M1 12h6M17 12h6M4.22 19.78l4.24-4.24M15.54 8.46l4.24-4.24"/>
                                </svg>
                                <span>Apple Pay, Google Pay & PayPal</span>
                            </div>
                            <div className="feature">
                                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                                    <rect x="2" y="4" width="20" height="16" rx="2"/>
                                    <path d="M7 15h0M17 15h0M7 11h0M13 11h0"/>
                                </svg>
                                <span>One-time setup, lifetime access</span>
                            </div>
                        </div>
                    </div>
                )}
            </div>
        </div>
    );
};

export default QuickSetup;
