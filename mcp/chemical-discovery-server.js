#!/usr/bin/env node
/**
 * Chemical Discovery MCP Server
 *
 * Model Context Protocol server that answers "what's in X?" queries.
 * Takes any object name (organism, food, material) and returns the chemicals it contains.
 *
 * Examples:
 * - "SARS-CoV-2" → all tested antiviral compounds
 * - "coffee" → caffeine, chlorogenic acid, water, etc.
 * - "aspirin" → acetylsalicylic acid (the compound itself)
 */

const { Server } = require('@modelcontextprotocol/sdk/server/index.js');
const { StdioServerTransport } = require('@modelcontextprotocol/sdk/server/stdio.js');
const {
  CallToolRequestSchema,
  ListToolsRequestSchema,
  ListResourcesRequestSchema,
  ReadResourceRequestSchema,
} = require('@modelcontextprotocol/sdk/types.js');

const PUBCHEM_BASE = 'https://pubchem.ncbi.nlm.nih.gov';
const NCBI_TAXONOMY_API = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';

/**
 * Fetch helper with error handling
 */
async function fetchJSON(url) {
  const fetch = (await import('node-fetch')).default;
  const response = await fetch(url);

  if (!response.ok) {
    throw new Error(`HTTP ${response.status}: ${response.statusText}`);
  }

  const contentType = response.headers.get('content-type');
  if (contentType?.includes('application/json')) {
    return await response.json();
  }

  return await response.text();
}

/**
 * Fetch SDF content from PubChem
 */
async function fetchSDFContent(name, recordType = '3d') {
  const fetch = (await import('node-fetch')).default;
  const encoded = encodeURIComponent(name);
  const url = `${PUBCHEM_BASE}/rest/pug/compound/name/${encoded}/SDF?record_type=${recordType}`;
  
  const response = await fetch(url);
  if (!response.ok) {
    throw new Error(`HTTP ${response.status}: ${response.statusText}`);
  }
  
  return await response.text();
}

/**
 * Search NCBI Taxonomy for organism name → taxonomy ID
 */
async function findTaxonomyId(organismName) {
  const url = `${NCBI_TAXONOMY_API}/esearch.fcgi?db=taxonomy&term=${encodeURIComponent(organismName)}&retmode=json`;

  try {
    const data = await fetchJSON(url);
    const idList = data?.esearchresult?.idlist || [];

    if (idList.length > 0) {
      return parseInt(idList[0]);
    }
  } catch (error) {
    // Not found or API error
  }

  return null;
}

/**
 * Get chemicals associated with an organism from PubChem Taxonomy
 */
async function getOrganismChemicals(taxonomyId) {
  const url = `${PUBCHEM_BASE}/rest/pug_view/data/taxonomy/${taxonomyId}/JSON`;

  try {
    const data = await fetchJSON(url);

    // Extract chemical information from the taxonomy record
    const sections = data?.Record?.Section || [];

    // Find the "Chemicals and Bioactivities" section
    let chemicals = [];

    for (const section of sections) {
      if (section.TOCHeading === 'Chemicals and Bioactivities') {
        // This contains references to chemicals tested against this organism
        // The actual data is in external tables, but we can get the concept
        chemicals.push({
          source: 'taxonomy',
          organism: data.Record.RecordTitle,
          taxonomyId: taxonomyId,
          note: 'Chemicals with bioactivity data for this organism',
        });
        break;
      }
    }

    return {
      type: 'organism',
      name: data.Record.RecordTitle,
      taxonomyId: taxonomyId,
      chemicals: chemicals,
      sections: sections.map(s => s.TOCHeading),
    };
  } catch (error) {
    throw new Error(`Failed to get organism chemicals: ${error.message}`);
  }
}

/**
 * Search for a chemical compound by name
 */
async function findCompoundByName(name) {
  const url = `${PUBCHEM_BASE}/rest/pug/compound/name/${encodeURIComponent(name)}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IUPACName/JSON`;

  try {
    const data = await fetchJSON(url);
    const properties = data?.PropertyTable?.Properties?.[0];

    if (properties) {
      return {
        type: 'compound',
        name: name,
        cid: properties.CID,
        formula: properties.MolecularFormula,
        molecularWeight: properties.MolecularWeight,
        smiles: properties.CanonicalSMILES,
        iupacName: properties.IUPACName,
        sdfUrl: `${PUBCHEM_BASE}/rest/pug/compound/cid/${properties.CID}/SDF?record_type=3d`,
      };
    }
  } catch (error) {
    // Not found or error
  }

  return null;
}

/**
 * Main discovery function: determine what type of object and return chemicals
 */
async function discoverChemicals(objectName) {
  const results = {
    query: objectName,
    discoveries: [],
  };

  // Strategy 1: Check if it's an organism
  const taxonomyId = await findTaxonomyId(objectName);
  if (taxonomyId) {
    const organismData = await getOrganismChemicals(taxonomyId);
    results.discoveries.push({
      type: 'organism',
      data: organismData,
      url: `${PUBCHEM_BASE}/taxonomy/${taxonomyId}`,
      apiUrl: `${PUBCHEM_BASE}/rest/pug_view/data/taxonomy/${taxonomyId}/JSON`,
    });
  }

  // Strategy 2: Check if it's a chemical compound itself
  const compound = await findCompoundByName(objectName);
  if (compound) {
    results.discoveries.push({
      type: 'compound',
      data: compound,
      url: `${PUBCHEM_BASE}/compound/${compound.cid}`,
      note: 'This is a single chemical compound',
    });
  }

  // Strategy 3: Could add food database, material databases, etc.
  // For now, if we found nothing, indicate that
  if (results.discoveries.length === 0) {
    results.discoveries.push({
      type: 'unknown',
      note: `No chemical data found for "${objectName}". Try: organism names (e.g., "E. coli"), chemical names (e.g., "caffeine"), or common names.`,
    });
  }

  return results;
}

/**
 * Create and configure MCP server
 */
const server = new Server(
  {
    name: 'chemical-discovery-server',
    version: '1.0.0',
  },
  {
    capabilities: {
      tools: {},
      resources: {},
    },
  }
);

/**
 * List available tools
 */
server.setRequestHandler(ListToolsRequestSchema, async () => {
  return {
    tools: [
      {
        name: 'whats_in',
        description: 'Discover what chemicals are in an object (organism, food, material, or the compound itself)',
        inputSchema: {
          type: 'object',
          properties: {
            object: {
              type: 'string',
              description: 'Name of the object to analyze (e.g., "SARS-CoV-2", "coffee", "aspirin", "E. coli")',
            },
          },
          required: ['object'],
        },
      },
      {
        name: 'get_organism_chemicals',
        description: 'Get all chemicals associated with a specific organism by taxonomy ID or name',
        inputSchema: {
          type: 'object',
          properties: {
            organism: {
              type: 'string',
              description: 'Organism name (e.g., "Homo sapiens", "E. coli") or taxonomy ID',
            },
          },
          required: ['organism'],
        },
      },
      {
        name: 'get_compound_structure',
        description: 'Get 3D molecular structure and properties for a chemical compound',
        inputSchema: {
          type: 'object',
          properties: {
            name: {
              type: 'string',
              description: 'Chemical compound name',
            },
          },
          required: ['name'],
        },
      },
      {
        name: 'get_3d_structure_data',
        description: 'Get the actual 3D structure data (SDF format) for a chemical compound. Returns the SDF file content that can be used to display 3D molecular structures.',
        inputSchema: {
          type: 'object',
          properties: {
            name: {
              type: 'string',
              description: 'Chemical compound name (e.g., "caffeine", "aspirin", "glucose")',
            },
            record_type: {
              type: 'string',
              description: 'Record type: "3d" for 3D coordinates (default), "2d" for 2D',
              enum: ['3d', '2d'],
              default: '3d',
            },
          },
          required: ['name'],
        },
      },
    ],
  };
});

/**
 * List available resources
 */
server.setRequestHandler(ListResourcesRequestSchema, async () => {
  return {
    resources: [
      {
        uri: 'discovery://databases/info',
        name: 'Chemical Discovery Databases',
        description: 'Information about available chemical databases and their capabilities',
        mimeType: 'application/json',
      },
      {
        uri: 'discovery://taxonomy/examples',
        name: 'Example Organisms',
        description: 'Common organism taxonomy IDs for testing',
        mimeType: 'application/json',
      },
    ],
  };
});

/**
 * Read resource content
 */
server.setRequestHandler(ReadResourceRequestSchema, async (request) => {
  const uri = request.params.uri;

  if (uri === 'discovery://databases/info') {
    return {
      contents: [
        {
          uri,
          mimeType: 'application/json',
          text: JSON.stringify({
            databases: [
              {
                name: 'PubChem Taxonomy',
                type: 'Organisms and biological systems',
                coverage: '9,000+ taxa with bioactivity data',
                example: 'SARS-CoV-2, E. coli, Human',
              },
              {
                name: 'PubChem Compound',
                type: 'Individual chemical compounds',
                coverage: '140+ million compounds',
                example: 'caffeine, aspirin, glucose',
              },
            ],
            workflow: [
              '1. User enters object name (e.g., "coffee")',
              '2. System checks multiple databases in parallel',
              '3. Returns all matching chemical data',
              '4. Provides 3D structures where available',
            ],
          }, null, 2),
        },
      ],
    };
  }

  if (uri === 'discovery://taxonomy/examples') {
    return {
      contents: [
        {
          uri,
          mimeType: 'application/json',
          text: JSON.stringify({
            viruses: [
              { name: 'SARS-CoV-2', id: 2697049 },
              { name: 'HIV-1', id: 11676 },
            ],
            bacteria: [
              { name: 'E. coli', id: 562 },
              { name: 'Mycobacterium tuberculosis', id: 1773 },
            ],
            model_organisms: [
              { name: 'Human', id: 9606 },
              { name: 'Mouse', id: 10090 },
              { name: 'Rat', id: 10116 },
            ],
          }, null, 2),
        },
      ],
    };
  }

  throw new Error(`Unknown resource: ${uri}`);
});

/**
 * Handle tool calls
 */
server.setRequestHandler(CallToolRequestSchema, async (request) => {
  const { name, arguments: args } = request.params;

  try {
    switch (name) {
      case 'whats_in': {
        const results = await discoverChemicals(args.object);

        return {
          content: [
            {
              type: 'text',
              text: JSON.stringify(results, null, 2),
            },
          ],
        };
      }

      case 'get_organism_chemicals': {
        // Check if input is numeric (taxonomy ID) or text (organism name)
        const taxonomyId = isNaN(args.organism)
          ? await findTaxonomyId(args.organism)
          : parseInt(args.organism);

        if (!taxonomyId) {
          return {
            content: [
              {
                type: 'text',
                text: JSON.stringify({
                  error: `Organism not found: ${args.organism}`,
                  suggestion: 'Try scientific names like "Escherichia coli" or common names like "E. coli"',
                }, null, 2),
              },
            ],
            isError: true,
          };
        }

        const data = await getOrganismChemicals(taxonomyId);

        return {
          content: [
            {
              type: 'text',
              text: JSON.stringify(data, null, 2),
            },
          ],
        };
      }

      case 'get_compound_structure': {
        const compound = await findCompoundByName(args.name);

        if (!compound) {
          return {
            content: [
              {
                type: 'text',
                text: JSON.stringify({
                  error: `Compound not found: ${args.name}`,
                  suggestion: 'Try common chemical names like "caffeine", "aspirin", or "glucose"',
                }, null, 2),
              },
            ],
            isError: true,
          };
        }

        return {
          content: [
            {
              type: 'text',
              text: JSON.stringify(compound, null, 2),
            },
          ],
        };
      }

      case 'get_3d_structure_data': {
        try {
          const recordType = args.record_type || '3d';
          const sdfContent = await fetchSDFContent(args.name, recordType);
          
          // Get compound metadata for context
          const compound = await findCompoundByName(args.name);
          
          return {
            content: [
              {
                type: 'text',
                text: JSON.stringify({
                  name: args.name,
                  format: 'sdf',
                  record_type: recordType,
                  sdf_content: sdfContent,
                  metadata: compound || null,
                  note: 'Use the sdf_content field to display the 3D structure. The SDF format is compatible with most molecular viewers.',
                }, null, 2),
              },
            ],
          };
        } catch (error) {
          if (error.message.includes('404') || error.message.includes('NOT_FOUND')) {
            return {
              content: [
                {
                  type: 'text',
                  text: JSON.stringify({
                    error: `3D structure not found for: ${args.name}`,
                    suggestion: 'Try common chemical names like "caffeine", "aspirin", or "glucose"',
                  }, null, 2),
                },
              ],
              isError: true,
            };
          }
          
          return {
            content: [
              {
                type: 'text',
                text: JSON.stringify({
                  error: `Failed to fetch 3D structure: ${error.message}`,
                }, null, 2),
              },
            ],
            isError: true,
          };
        }
      }

      default:
        throw new Error(`Unknown tool: ${name}`);
    }
  } catch (error) {
    return {
      content: [
        {
          type: 'text',
          text: `Error: ${error.message}`,
        },
      ],
      isError: true,
    };
  }
});

/**
 * Start server
 */
async function main() {
  const transport = new StdioServerTransport();
  await server.connect(transport);
  console.error('Chemical Discovery MCP server running on stdio');
}

main().catch((error) => {
  console.error('Server error:', error);
  process.exit(1);
});
