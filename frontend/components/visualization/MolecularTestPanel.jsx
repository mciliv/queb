import React, { useState } from 'react';

const MolecularTestPanel = ({ onTestAnalysis, isProcessing }) => {
  const [showTestPanel, setShowTestPanel] = useState(false);

  // Test objects using minimal subset validation approach
  const testObjects = [
    // Basics - high confidence minimal subset
    { name: '💧 Water', input: 'water', expected: ['water'], category: 'basics' },
    { name: '🧪 Ethanol', input: 'ethanol', expected: ['ethanol'], category: 'basics' },
    { name: '🧂 Salt', input: 'sodium chloride', expected: ['sodium', 'chloride'], category: 'basics' },
    
    // Beverages - essential components only (minimal subset)
    { name: '🍷 Wine', input: 'red wine', expected: ['ethanol', 'water'], category: 'beverages' },
    { name: '☕ Coffee', input: 'black coffee', expected: ['water', 'caffeine'], category: 'beverages' },
    
    // Biological - most confident components only (minimal subset)  
    { name: '🍎 Apple', input: 'fresh apple', expected: ['water'], category: 'biological' },
    { name: '🥬 Kale', input: 'kale', expected: ['water'], category: 'biological' },
    { name: '🥛 Milk', input: 'milk', expected: ['water'], category: 'biological' }
  ];

  const handleTestClick = async (testObject) => {
    if (isProcessing) return;
    
    console.log(`🧪 Testing molecular analysis for: ${testObject.name}`);
    console.log(`   Input: "${testObject.input}"`);
    console.log(`   Category: ${testObject.category}`);
    console.log(`   Expected subset: ${testObject.expected.join(', ')} (minimal required components)`);
    
    await onTestAnalysis(testObject.input);
  };

  const runAllTests = async () => {
    if (isProcessing) return;
    
    console.log('🧪 Running all molecular visualization tests in main interface...');
    console.log('   Each test will inject into the main UI for visual validation');
    
    for (const testObject of testObjects) {
      console.log(`\n🔬 Testing: ${testObject.name}`);
      console.log(`   Category: ${testObject.category} | Expected: ${testObject.expected.join(', ')}`);
      
      await handleTestClick(testObject);
      
      // Wait longer for visual inspection of results
      await new Promise(resolve => setTimeout(resolve, 4000)); 
    }
    
    console.log('✅ All visual molecular tests completed - check main interface results');
  };

  if (!showTestPanel) {
    return (
      <button 
        className="test-panel-toggle"
        onClick={() => setShowTestPanel(true)}
        title="Open molecular test panel"
        style={{
          position: 'fixed',
          bottom: '20px',
          left: '20px',
          background: 'rgba(0, 0, 0, 0.7)',
          color: 'white',
          border: 'none',
          borderRadius: '50%',
          width: '50px',
          height: '50px',
          fontSize: '20px',
          cursor: 'pointer',
          zIndex: 1000
        }}
      >
        🧪
      </button>
    );
  }

  return (
    <div 
      className="molecular-test-panel"
      style={{
        position: 'fixed',
        bottom: '20px',
        left: '20px',
        background: 'rgba(0, 0, 0, 0.9)',
        color: 'white',
        padding: '20px',
        borderRadius: '10px',
        maxWidth: '400px',
        maxHeight: '70vh',
        overflow: 'auto',
        zIndex: 1000
      }}
    >
      <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '15px' }}>
        <h3 style={{ margin: 0 }}>🧪 Molecular Test Panel</h3>
        <button 
          onClick={() => setShowTestPanel(false)}
          style={{ background: 'none', border: 'none', color: 'white', fontSize: '20px', cursor: 'pointer' }}
        >
          ×
        </button>
      </div>

      <div style={{ marginBottom: '15px' }}>
        <button
          onClick={runAllTests}
          disabled={isProcessing}
          style={{
            background: '#4CAF50',
            color: 'white',
            border: 'none',
            padding: '10px 15px',
            borderRadius: '5px',
            cursor: isProcessing ? 'not-allowed' : 'pointer',
            width: '100%',
            marginBottom: '10px'
          }}
        >
          {isProcessing ? 'Running Tests...' : '🚀 Run All Tests'}
        </button>
      </div>

      <div style={{ marginBottom: '15px' }}>
        <h4 style={{ margin: '0 0 10px 0', fontSize: '14px' }}>Quick Tests:</h4>
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(120px, 1fr))', gap: '8px' }}>
          {testObjects.map((testObject, index) => (
            <button
              key={index}
              onClick={() => handleTestClick(testObject)}
              disabled={isProcessing}
              style={{
                background: getCategoryColor(testObject.category),
                color: 'white',
                border: 'none',
                padding: '8px 12px',
                borderRadius: '5px',
                cursor: isProcessing ? 'not-allowed' : 'pointer',
                fontSize: '12px',
                textAlign: 'center'
              }}
              title={`Expected subset: ${testObject.expected.join(', ')}`}
            >
              {testObject.name}
            </button>
          ))}
        </div>
      </div>

      <div style={{ fontSize: '12px', opacity: 0.7 }}>
        <p style={{ margin: '5px 0' }}>💡 Click to inject test case into main interface</p>
        <p style={{ margin: '5px 0' }}>🏃 Run All Tests to see visual comparison across all cases</p>  
        <p style={{ margin: '5px 0' }}>🎯 Uses minimal subset validation - additional compounds are expected</p>
      </div>
    </div>
  );
};

const getCategoryColor = (category) => {
  const colors = {
    basics: '#2196F3',      // Blue - high confidence 
    beverages: '#FF9800',   // Orange - moderate confidence
    biological: '#4CAF50'   // Green - basic confidence
  };
  return colors[category] || '#757575';
};

export default MolecularTestPanel;