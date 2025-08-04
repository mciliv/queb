// Manual mock for payment.js
export const paymentManager = {
  checkPaymentMethod: jest.fn().mockResolvedValue(true),
  incrementUsage: jest.fn().mockResolvedValue()
};