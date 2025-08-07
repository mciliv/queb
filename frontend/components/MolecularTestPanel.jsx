import React, { useState } from 'react';

const MolecularTestPanel = ({ onTestAnalysis, isProcessing }) => {
  const [showTestPanel, setShowTestPanel] = useState(false);

  // Predefined test objects for quick molecular analysis testing
  const testObjects = [
    {
      name: "Water",
      description: "water",
      expectedMolecules: ["H2O"],
      category: "simple"
    },
    {
      name: "Apple",
      description: "apple",
      expectedMolecules: ["water", "glucose", "fructose", "cellulose"],
      category: "food"
    },
    {
      name: "Coffee",
      description: "coffee",
      expectedMolecules: ["caffeine", "water", "acids"],
      category: "beverage"
    },
    {
      name: "Kale",
      description: "kale",
      expectedMolecules: ["water", "glucose", "chlorophyll", "cellulose"],
      category: "vegetable"
    },
    {
      name: "Ethanol",
      description: "ethanol",
      expectedMolecules: ["C2H6O"],
      category: "chemical"
    },
    {
      name: "Salt",
      description: "table salt",
      expectedMolecules: ["NaCl"],
      category: "simple"
    },
    {
      name: "Wine",
      description: "red wine",
      expectedMolecules: ["ethanol", "water", "tartaric acid", "resveratrol"],
      category: "beverage"
    },
    {
      name: "Plastic Bottle",
      description: "plastic water bottle",
      expectedMolecules: ["polyethylene", "additives"],
      category: "material"
    }
  ];

  const handleTestClick = async (testObject) => {
    if (isProcessing) return;
    
    console.log(`🧪 Testing molecular analysis for: ${testObject.name}`);
    console.log(`Expected molecules: ${testObject.expectedMolecules.join(', ')}`);
    
    await onTestAnalysis(testObject.description);
  };

  const runAllTests = async () => {
    if (isProcessing) return;
    
    console.log('🧪 Running all molecular visualization tests...');
    
    for (const testObject of testObjects) {
      console.log(`\n🔬 Testing: ${testObject.name}`);
      await new Promise(resolve => setTimeout(resolve, 1000)); // Brief delay between tests
      await handleTestClick(testObject);
      await new Promise(resolve => setTimeout(resolve, 3000)); // Wait for analysis to complete
    }
    
    console.log('✅ All molecular tests completed');
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
              title={`Expected: ${testObject.expectedMolecules.join(', ')}`}
            >
              {testObject.name}
            </button>
          ))}
        </div>
      </div>

      <div style={{ fontSize: '12px', opacity: 0.7 }}>
        <p style={{ margin: '5px 0' }}>💡 Click any test to analyze that object</p>
        <p style={{ margin: '5px 0' }}>🏃 Run All Tests for comprehensive testing</p>
        <p style={{ margin: '5px 0' }}>🎯 Check console for expected vs actual results</p>
      </div>
    </div>
  );
};

const getCategoryColor = (category) => {
  const colors = {
    simple: '#2196F3',    // Blue
    food: '#4CAF50',      // Green  
    beverage: '#FF9800',  // Orange
    chemical: '#9C27B0',  // Purple
    material: '#607D8B',  // Blue Grey
    vegetable: '#8BC34A'  // Light Green
  };
  return colors[category] || '#757575';
};

export default MolecularTestPanel;