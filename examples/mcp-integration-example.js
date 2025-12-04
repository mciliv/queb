/**
 * MCP Integration Examples
 * 
 * Shows how to use MCP-enhanced services in your existing code
 */

const PerformanceCache = require('../src/server/services/performance-cache');
const { getMCPService } = require('../src/server/services/mcp-integration');

// ============================================================================
// EXAMPLE 1: Add performance tracking and caching to text analysis
// ============================================================================

async function analyzeTextWithCache(description) {
  const cache = new PerformanceCache({ logger: console });
  
  // Wrap expensive operation
  const result = await cache.wrap(
    'analyzeText',          // Operation name
    description,            // Cache key input
    async () => {
      // Your expensive operation here
      console.log('Calling AI service (expensive)...');
      
      // Simulate AI call
      await new Promise(resolve => setTimeout(resolve, 2000));
      
      return {
        object: description,
        molecules: [
          { name: 'Water', smiles: 'O', formula: 'H2O' }
        ]
      };
    },
    {
      ttl: 7 * 24 * 60 * 60 * 1000,  // 7 days
      metadata: {
        inputLength: description.length,
        source: 'user_input'
      }
    }
  );
  
  return result;
}

// Usage
async function example1() {
  console.log('\n=== Example 1: Cached Text Analysis ===\n');
  
  // First call: cache miss, slow
  console.log('First call:');
  await analyzeTextWithCache('coffee');
  
  // Second call: cache hit, instant
  console.log('\nSecond call (should be cached):');
  await analyzeTextWithCache('coffee');
  
  // Get statistics
  const cache = new PerformanceCache();
  const stats = await cache.getStats('analyzeText');
  console.log('\nStatistics:', stats);
}

// ============================================================================
// EXAMPLE 2: Use MCP service for file operations
// ============================================================================

async function listMolecularFiles() {
  const mcpService = getMCPService({ logger: console });
  
  console.log('\n=== Example 2: List SDF Files ===\n');
  
  // List all SDF files with metadata
  const sdfFiles = await mcpService.listSDFFiles('src/backend/sdf_files');
  
  console.log(`Found ${sdfFiles.length} SDF files:\n`);
  
  sdfFiles.forEach(file => {
    console.log(`  📄 ${file.name}`);
    console.log(`     Size: ${(file.size / 1024).toFixed(2)} KB`);
    console.log(`     Modified: ${file.modified.toLocaleString()}`);
    console.log('');
  });
  
  return sdfFiles;
}

// ============================================================================
// EXAMPLE 3: Database schema inspection
// ============================================================================

async function inspectDatabase() {
  const mcpService = getMCPService({ logger: console });
  
  console.log('\n=== Example 3: Database Inspection ===\n');
  
  try {
    // Get all tables
    const tables = await mcpService.getDatabaseSchema();
    
    console.log(`Database has ${tables.length} tables:\n`);
    
    for (const table of tables) {
      console.log(`📊 ${table.table_name} (${table.table_type})`);
      
      // Get columns for this table
      const columns = await mcpService.getDatabaseSchema(table.table_name);
      columns.forEach(col => {
        const nullable = col.is_nullable === 'YES' ? 'NULL' : 'NOT NULL';
        console.log(`   - ${col.column_name}: ${col.data_type} ${nullable}`);
      });
      console.log('');
    }
  } catch (error) {
    console.log('Database not configured or not accessible');
    console.log('This is optional - the app works without a database');
  }
}

// ============================================================================
// EXAMPLE 4: Performance monitoring
// ============================================================================

async function monitorPerformance() {
  const cache = new PerformanceCache();
  
  console.log('\n=== Example 4: Performance Monitoring ===\n');
  
  // Simulate multiple operations
  const operations = ['analyzeText', 'analyzeImage', 'generateSDF'];
  
  for (const op of operations) {
    const stats = await cache.getStats(op, 24);
    
    if (stats) {
      console.log(`${op}:`);
      console.log(`  Total calls: ${stats.totalCalls}`);
      console.log(`  Cache hit rate: ${stats.cacheHitRate}`);
      console.log(`  Avg duration: ${stats.avgDuration}`);
      console.log(`  Min/Max: ${stats.minDuration} / ${stats.maxDuration}`);
      console.log('');
    } else {
      console.log(`${op}: No data yet\n`);
    }
  }
}

// ============================================================================
// EXAMPLE 5: Integrate into existing service
// ============================================================================

class EnhancedMolecularService {
  constructor() {
    this.cache = new PerformanceCache({ 
      logger: console,
      enableCache: true,
      trackPerformance: true
    });
    this.mcpService = getMCPService();
  }

  async analyzeMolecule(input) {
    // Automatic caching and performance tracking
    return await this.cache.wrap(
      'analyzeMolecule',
      input,
      async () => {
        // Your existing analysis logic
        return this._performAnalysis(input);
      }
    );
  }

  async _performAnalysis(input) {
    // Simulate expensive operation
    console.log(`Analyzing: ${input}`);
    await new Promise(resolve => setTimeout(resolve, 1000));
    return { result: 'analysis complete' };
  }

  async getPerformanceReport() {
    return await this.cache.getStats('analyzeMolecule', 24);
  }

  async getCacheEfficiency() {
    return await this.cache.getCacheStats();
  }
}

// ============================================================================
// Run examples
// ============================================================================

async function runAllExamples() {
  console.log('\n🧬 MCP Integration Examples\n');
  console.log('These examples show how to enhance your existing code with MCP services.\n');
  
  try {
    // Example 1: Caching
    await example1();
    
    // Example 2: File operations
    await listMolecularFiles();
    
    // Example 3: Database inspection
    await inspectDatabase();
    
    // Example 4: Performance monitoring
    await monitorPerformance();
    
    // Example 5: Enhanced service
    console.log('\n=== Example 5: Enhanced Service ===\n');
    const service = new EnhancedMolecularService();
    await service.analyzeMolecule('caffeine');
    await service.analyzeMolecule('caffeine'); // Should be cached
    
    const report = await service.getPerformanceReport();
    console.log('\nPerformance Report:', report);
    
    console.log('\n✅ All examples completed!\n');
    
  } catch (error) {
    console.error('❌ Error running examples:', error.message);
    console.log('\nNote: Some examples require database setup (optional)');
    console.log('Run: node database/setup/setup-database.js');
    console.log('Then: psql -U mol_user -d mol_users -f database/migrations/001_create_mcp_tables.sql\n');
  }
  
  process.exit(0);
}

// Run if called directly
if (require.main === module) {
  runAllExamples();
}

module.exports = {
  analyzeTextWithCache,
  listMolecularFiles,
  inspectDatabase,
  monitorPerformance,
  EnhancedMolecularService
};









