import React from 'react';
import { useStripeSetup } from '../hooks/useStripeSetup';

const AccountSetup = () => {
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