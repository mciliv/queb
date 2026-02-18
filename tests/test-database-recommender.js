#!/usr/bin/env node

/**
 * Test script for database recommender and chemical contents finder
 * 
 * Usage: node tests/test-database-recommender.js [item-name] [mode]
 * Example: node tests/test-database-recommender.js "coffee" "query"
 * 
 * Modes:
 *   "recommend" (default) - Get database recommendations only
 *   "query" - Use AI to pick database and perform actual query
 * 
 * Note: Uses Node.js built-in fetch (available in Node 18+)
 */

const BASE_URL = process.env.API_URL || 'http://localhost:8080';
const ITEM = process.argv[2] || 'coffee';
const MODE = process.argv[3] || 'recommend';

async function testDatabaseRecommender() {
  console.log(`🧪 Testing Database ${MODE === 'query' ? 'Query' : 'Recommender'}\n`);
  console.log(`Item: "${ITEM}"`);
  console.log(`Mode: ${MODE}\n`);

  const endpoint = MODE === 'query' 
    ? `${BASE_URL}/api/find-chemical-contents`
    : `${BASE_URL}/api/recommend-database`;

  try {
    const response = await fetch(endpoint, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ item: ITEM }),
    });

    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`HTTP ${response.status}: ${errorText}`);
    }

    const result = await response.json();

    if (MODE === 'query') {
      // Display query results
      console.log('✅ Success!\n');
      console.log('Selected Database:');
      console.log(`  Name: ${result.selectedDatabase.name}`);
      console.log(`  Reason: ${result.selectedDatabase.reason}`);
      console.log(`  Query URL: ${result.queryInfo.url}\n`);

      if (result.chemicals && result.chemicals.length > 0) {
        console.log(`Found ${result.chemicals.length} chemical(s):\n`);
        result.chemicals.forEach((chem, index) => {
          console.log(`  ${index + 1}. ${chem.name || 'Unknown'}`);
          if (chem.formula) console.log(`     Formula: ${chem.formula}`);
          if (chem.molecularWeight) console.log(`     Molecular Weight: ${chem.molecularWeight}`);
          if (chem.smiles) console.log(`     SMILES: ${chem.smiles}`);
          if (chem.iupac) console.log(`     IUPAC: ${chem.iupac}`);
          if (chem.cid) console.log(`     PubChem CID: ${chem.cid}`);
          console.log('');
        });
      } else {
        console.log('⚠️  No chemicals found for this item.\n');
      }

      if (result.metadata) {
        console.log('Metadata:');
        console.log(`  Source: ${result.metadata.source}`);
        console.log(`  Query Type: ${result.metadata.queryType}`);
        console.log(`  Result Count: ${result.metadata.resultCount || 0}\n`);
      }
    } else {
      // Display recommendations
      console.log('✅ Success!\n');
      console.log('Analysis:');
      console.log(`  Type: ${result.analysis.type}`);
      console.log(`  Confidence: ${result.analysis.confidence}`);
      console.log(`  Characteristics: ${result.analysis.characteristics.join(', ')}\n`);

      console.log('Recommendations:');
      result.recommendations.forEach((rec, index) => {
        console.log(`  ${index + 1}. ${rec.name} (Priority: ${rec.priority})`);
        console.log(`     Reason: ${rec.reason}`);
      });

      console.log('\nDatabase Details:');
      result.databases.forEach((db, index) => {
        console.log(`\n  ${index + 1}. ${db.name}`);
        console.log(`     Provider: ${db.provider}`);
        console.log(`     Description: ${db.description}`);
        console.log(`     Base URL: ${db.baseUrl}`);
        console.log(`     Authentication: ${db.authentication}`);
        if (db.rateLimit) {
          console.log(`     Rate Limit: ${db.rateLimit}`);
        }
      });
      console.log('\n');
    }
  } catch (error) {
    console.error('❌ Error:', error.message);
    if (error.cause) {
      console.error('   Cause:', error.cause);
    }
    
    // Provide helpful error messages
    if (error.message.includes('ECONNREFUSED') || error.message.includes('fetch failed')) {
      console.error('\n💡 Tip: Make sure the server is running:');
      console.error('   npm start');
      console.error('   or');
      console.error('   npm run debug');
    }
    
    process.exit(1);
  }
}

testDatabaseRecommender();
