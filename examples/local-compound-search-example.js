/**
 * local-compound-search-example.js
 *
 * Example usage of the LocalCompoundSearch service for searching compounds
 * based on object names with fuzzy matching and local caching.
 */

const LocalCompoundSearch = require('../src/server/services/local-compound-search');

async function main() {
  console.log('=== Local Compound Search Example ===\n');

// Initialize the search service
const searchService = new LocalCompoundSearch({
  maxResults: 10,
  minScore: 0.3,  // Lower threshold for more results
  fuzzyThreshold: 0.6  // More lenient fuzzy matching
});

  console.log('1. Initial database stats:');
  console.log(searchService.getStats());
  console.log();

  // Example 1: Add some compounds manually
  console.log('2. Adding sample compounds...');
  const sampleCompounds = [
    {
      cid: '2244',
      name: 'Aspirin',
      names: ['Acetylsalicylic acid', 'ASA', '2-Acetoxybenzoic acid'],
      smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O',
      iupac: '2-acetyloxybenzoic acid',
      formula: 'C9H8O4',
      molecularWeight: 180.16,
      source: 'manual'
    },
    {
      cid: '2519',
      name: 'Caffeine',
      names: ['1,3,7-Trimethylpurine-2,6-dione', 'Guaranine'],
      smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
      iupac: '1,3,7-trimethylpurine-2,6-dione',
      formula: 'C8H10N4O2',
      molecularWeight: 194.19,
      source: 'manual'
    },
    {
      cid: '702',
      name: 'Ethanol',
      names: ['Alcohol', 'Ethyl alcohol', 'Grain alcohol'],
      smiles: 'CCO',
      iupac: 'ethanol',
      formula: 'C2H6O',
      molecularWeight: 46.07,
      source: 'manual'
    }
  ];

  sampleCompounds.forEach(compound => {
    searchService.addCompound(compound);
  });
  console.log(`Added ${sampleCompounds.length} compounds to local database\n`);

  // Example 2: Search for compounds
  console.log('3. Searching for compounds...');

  const searchQueries = [
    'aspirin',      // Exact match
    'caffine',      // Fuzzy match (typo)
    'alcohol',      // Alternative name
    'acetic acid',  // Partial match
    'C8H10N4O2'     // Formula search
  ];

  searchQueries.forEach(query => {
    console.log(`\nSearch for "${query}":`);
    const results = searchService.search(query);

    if (results.length === 0) {
      console.log('  No results found');
    } else {
      results.forEach((result, index) => {
        console.log(`  ${index + 1}. ${result.compound.name} (CID: ${result.compound.cid})`);
        console.log(`     Score: ${(result.score * 100).toFixed(1)}%, Match: ${result.matchType}`);
        console.log(`     Formula: ${result.compound.formula}, MW: ${result.compound.molecularWeight}`);
      });
    }
  });

  // Example 3: Get compound by CID
  console.log('\n4. Getting compound by CID...');
  const aspirin = searchService.getCompoundByCID('2244');
  if (aspirin) {
    console.log(`Found: ${aspirin.name}`);
    console.log(`SMILES: ${aspirin.smiles}`);
    console.log(`IUPAC: ${aspirin.iupac}`);
  }

  // Example 4: Populate from PubChem (if network available)
  console.log('\n5. Populating from PubChem...');
  try {
    const pubchemCompounds = ['glucose', 'ibuprofen', 'paracetamol', 'penicillin'];
    console.log(`Fetching ${pubchemCompounds.length} compounds from PubChem...`);

    const populateResults = await searchService.populateFromPubChem(pubchemCompounds);
    console.log('Population results:', populateResults);

    if (populateResults.success > 0) {
      // Search again to show new results
      console.log('\nSearch for "glucose" after PubChem population:');
      const glucoseResults = searchService.search('glucose');
      glucoseResults.forEach((result, index) => {
        console.log(`  ${index + 1}. ${result.compound.name} (CID: ${result.compound.cid})`);
      });
    } else {
      console.log('No compounds were successfully added from PubChem');
    }

  } catch (error) {
    console.log('PubChem population failed:', error.message);
    console.log('Continuing with local compounds only...');
  }

  // Example 5: Save and reload data
  console.log('\n6. Saving data and testing persistence...');
  searchService.saveLocalData();

  // Create a new instance to test loading
  const newSearchService = new LocalCompoundSearch();
  console.log('Reloaded database stats:', newSearchService.getStats());

  const reloadTest = newSearchService.search('aspirin');
  console.log(`Search for "aspirin" after reload: ${reloadTest.length} results`);

  // Example 6: Advanced search options
  console.log('\n7. Advanced search options...');

  // Search with different thresholds
  const strictResults = searchService.search('asp', { minScore: 0.8 });
  const lenientResults = searchService.search('asp', { minScore: 0.3, maxResults: 5 });

  console.log(`Strict search (minScore: 0.8): ${strictResults.length} results`);
  console.log(`Lenient search (minScore: 0.3): ${lenientResults.length} results`);

  console.log('\n=== Example Complete ===');

  // Cleanup example (commented out - remove to actually clear data)
  // console.log('\nClearing local database...');
  // searchService.clear();
}

// Run the example
if (require.main === module) {
  main().catch(console.error);
}

module.exports = { main };
