# Local Compound Search Service

A local search service for chemical compounds that provides fuzzy matching and fast retrieval without requiring external API calls.

## Features

- **Local Database**: Maintains a local cache of compound information
- **Fuzzy Search**: Uses Levenshtein distance for approximate string matching
- **Multiple Data Sources**: Can populate from PubChem, ChEMBL, and local files
- **Search Ranking**: Results are ranked by relevance score
- **Persistent Storage**: Saves data to JSON files for reuse across sessions
- **Performance**: Fast in-memory search with optional caching

## Quick Start

```javascript
const LocalCompoundSearch = require('./src/server/services/local-compound-search');

// Initialize the service
const search = new LocalCompoundSearch({
  maxResults: 10,      // Maximum results per search
  minScore: 0.4,       // Minimum similarity score (0-1)
  fuzzyThreshold: 0.7  // Fuzzy matching threshold
});

// Add some compounds
search.addCompound({
  cid: '2244',
  name: 'Aspirin',
  names: ['Acetylsalicylic acid', 'ASA'],
  smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O',
  formula: 'C9H8O4',
  molecularWeight: 180.16
});

// Search for compounds
const results = search.search('aspirin');
console.log(results[0].compound.name); // "Aspirin"
console.log(results[0].score); // 1.0 (exact match)
```

## API Reference

### Constructor Options

- `dataDir`: Directory for storing local data files (default: `data/`)
- `maxResults`: Maximum number of search results (default: 50)
- `minScore`: Minimum similarity score for results (default: 0.3)
- `fuzzyThreshold`: Threshold for fuzzy matching (default: 0.8)
- `enableCaching`: Whether to save/load data from files (default: true)

### Methods

#### `addCompound(compound)`

Adds a compound to the local database.

**Parameters:**
- `compound.cid` (required): PubChem Compound ID
- `compound.name`: Primary compound name
- `compound.names`: Array of alternative names
- `compound.smiles`: SMILES notation
- `compound.iupac`: IUPAC name
- `compound.formula`: Molecular formula
- `compound.molecularWeight`: Molecular weight

#### `search(query, options)`

Searches for compounds matching the query.

**Parameters:**
- `query`: Search string
- `options`: Search options (overrides constructor options)

**Returns:** Array of results with `compound`, `score`, and `matchType` properties

#### `getCompoundByCID(cid)`

Retrieves a compound by its PubChem CID.

#### `populateFromPubChem(compoundNames)`

Populates the local database by fetching compounds from PubChem.

**Parameters:**
- `compoundNames`: Array of compound names to fetch

#### `populateFromFile(filePath)`

Loads compounds from a JSON file.

#### `saveLocalData()`

Saves the current database to disk.

#### `getStats()`

Returns database statistics.

#### `clear()`

Clears all local data.

## Search Features

### Match Types

- **exact**: Exact string match (highest score)
- **fuzzy**: Approximate string match using Levenshtein distance
- **word**: Partial word match for multi-word queries

### Search Behavior

- Searches across primary names, IUPAC names, alternative names, and formulas
- Uses fuzzy matching for typos and variations
- Ranks results by relevance score
- Supports partial word matching

### Example Searches

```javascript
// Exact match
search.search('aspirin') // Score: 1.0

// Fuzzy match (typo)
search.search('caffine') // Matches "caffeine" with score ~0.9

// Alternative name
search.search('alcohol') // Matches "ethanol" via alternative names

// Formula search
search.search('C8H10N4O2') // Matches caffeine by formula

// Partial match
search.search('acetic') // Matches aspirin via "acetylsalicylic acid"
```

## Data Sources

### PubChem Integration

Automatically fetches compound data from PubChem using the existing name resolver:

```javascript
const results = await search.populateFromPubChem([
  'aspirin',
  'caffeine',
  'glucose',
  'ibuprofen'
]);
console.log(`Added ${results.success} compounds`);
```

### File Import

Load compounds from JSON files:

```javascript
search.populateFromFile('./my-compounds.json');
```

Expected JSON format:
```json
[
  {
    "cid": "2244",
    "name": "Aspirin",
    "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "formula": "C9H8O4"
  }
]
```

## Performance

- **Memory Usage**: Stores compounds in memory with optional file persistence
- **Search Speed**: Sub-millisecond searches for typical databases
- **Scalability**: Handles thousands of compounds efficiently
- **Caching**: Automatic persistence to avoid re-fetching data

## Integration with Queb

This service integrates with the existing Queb molecular processing pipeline:

```javascript
const { resolveName } = require('./name-resolver');
const MolecularProcessor = require('./molecular-processor');

// First, try local search
const localResults = search.search('aspirin');
if (localResults.length > 0) {
  const compound = localResults[0].compound;
  console.log(`Found locally: ${compound.name}`);
} else {
  // Fall back to PubChem
  const resolved = await resolveName('aspirin');
  search.addCompound(resolved);
}
```

## Running the Example

```bash
node examples/local-compound-search-example.js
```

This will demonstrate:
- Adding compounds manually
- Searching with different query types
- Populating from PubChem
- Data persistence
- Advanced search options
