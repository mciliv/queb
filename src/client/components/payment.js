/**
 * Payment Manager - Simple implementation for testing
 */

const paymentManager = {
  processPayment: async (amount) => {
    return { success: true, transactionId: `tx_${Date.now()}`, amount };
  },
  
  validateCard: (cardInfo) => {
    return { valid: true, errors: [] };
  },
  
  getPaymentMethods: () => {
    return ['card', 'paypal'];
  }
};

module.exports = {
  paymentManager
};


























