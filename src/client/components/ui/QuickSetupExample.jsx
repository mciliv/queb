import React from 'react';
import QuickSetup from './QuickSetup';

/**
 * Example usage of the QuickSetup component
 * This shows how to integrate the converted React component
 */
const QuickSetupExample = () => {
  return (
    <div>
      {/* Load Stripe script - required for payment functionality */}
      <script src="https://js.stripe.com/v3/"></script>
      
      {/* The converted QuickSetup component */}
      <QuickSetup />
    </div>
  );
};

export default QuickSetupExample;
