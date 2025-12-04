/**
 * Unit tests for MolecularProcessor
 * 
 * This test demonstrates the new test structure and patterns:
 * - Using test builders for consistent test data
 * - Custom matchers for domain-specific assertions  
 * - Proper mocking and isolation
 * - Clear test descriptions
 */

const MolecularProcessor = require('@/server/services/molecular-processor');
const { create } = require('@test/builders');
const { sleep, validators } = require('@test/helpers');
const { extendExpect } = require('@test/matchers');

// Extend Jest with custom matchers
extendExpect();

// Mock external dependencies
jest.mock('openai');
jest.mock('fs');

describe('MolecularProcessor', () => {
  let processor;
  let mockOpenAI;
  
  beforeEach(() => {
    // Reset all mocks
    jest.clearAllMocks();
    
    // Create fresh instance
    processor = new MolecularProcessor({
      apiKey: 'test-key',
      model: 'gpt-4'
    });
    
    // Get mock instance
    mockOpenAI = require('openai').OpenAI;
  });
  
  describe('constructor', () => {
    it('should initialize with required configuration', () => {
      expect(processor).toBeDefined();
      expect(processor.apiKey).toBe('test-key');
      expect(processor.model).toBe('gpt-4');
    });
    
    it('should throw error when API key is missing', () => {
      expect(() => {
        new MolecularProcessor({});
      }).toThrow('API key is required');
    });
  });
  
  describe('processText', () => {
    it('should extract molecules from text successfully', async () => {
      // Arrange
      const mockResponse = create.response()
        .withData({
          object: 'coffee',
          chemicals: [
            create.molecule().withName('Caffeine').withSmiles('CN1C=NC2=C1C(=O)N(C(=O)N2C)C').build(),
            create.molecule().withName('Chlorogenic acid').withSmiles('O=C(O)C1=CC(O)=C(O)C=C1').build()
          ]
        })
        .build();
      
      mockOpenAI.mockImplementation(() => ({
        chat: {
          completions: {
            create: jest.fn().mockResolvedValue({
              choices: [{
                message: {
                  content: JSON.stringify(mockResponse.data)
                }
              }]
            })
          }
        }
      }));
      
      // Act
      const result = await processor.processText('What chemicals are in coffee?');
      
      // Assert
      expect(result).toBeSuccessResponse();
      expect(result.data.chemicals).toContainMolecules(2);
      expect(result.data.object).toBe('coffee');
      
      // Verify all molecules have valid SMILES
      result.data.chemicals.forEach(molecule => {
        expect(molecule.smiles).toBeValidSmiles();
      });
    });
    
    it('should handle empty results gracefully', async () => {
      // Arrange
      mockOpenAI.mockImplementation(() => ({
        chat: {
          completions: {
            create: jest.fn().mockResolvedValue({
              choices: [{
                message: {
                  content: JSON.stringify({
                    object: 'unknown',
                    chemicals: []
                  })
                }
              }]
            })
          }
        }
      }));
      
      // Act
      const result = await processor.processText('Random text with no chemicals');
      
      // Assert
      expect(result).toBeSuccessResponse();
      expect(result.data.chemicals).toEqual([]);
    });
    
    it('should handle API errors appropriately', async () => {
      // Arrange
      mockOpenAI.mockImplementation(() => ({
        chat: {
          completions: {
            create: jest.fn().mockRejectedValue(new Error('API rate limit exceeded'))
          }
        }
      }));
      
      // Act
      const result = await processor.processText('Some text');
      
      // Assert
      expect(result).toBeErrorResponse(500);
      expect(result.message).toContain('rate limit');
    });
    
    it('should validate input text', async () => {
      // Act & Assert
      await expect(processor.processText('')).rejects.toThrow('Text cannot be empty');
      await expect(processor.processText(null)).rejects.toThrow('Text must be a string');
      await expect(processor.processText(123)).rejects.toThrow('Text must be a string');
    });
    
    it('should complete within performance threshold', async () => {
      // Arrange
      mockOpenAI.mockImplementation(() => ({
        chat: {
          completions: {
            create: jest.fn().mockImplementation(async () => {
              await sleep(100); // Simulate API latency
              return {
                choices: [{
                  message: {
                    content: JSON.stringify({
                      object: 'test',
                      chemicals: [create.molecule().build()]
                    })
                  }
                }]
              };
            })
          }
        }
      }));
      
      // Act & Assert
      await expect(async () => {
        await processor.processText('Quick test');
      }).toCompleteWithin(500);
    });
  });
  
  describe('processImage', () => {
    it('should process image and extract molecules', async () => {
      // Arrange
      const imageFile = create.file()
        .asPng()
        .withSize(1024 * 500) // 500KB
        .build();
      
      const mockResponse = create.response()
        .withData({
          object: 'chemical structure',
          chemicals: [
            create.molecule().asWater().build()
          ]
        })
        .build();
      
      mockOpenAI.mockImplementation(() => ({
        chat: {
          completions: {
            create: jest.fn().mockResolvedValue({
              choices: [{
                message: {
                  content: JSON.stringify(mockResponse.data)
                }
              }]
            })
          }
        }
      }));
      
      // Act
      const result = await processor.processImage(imageFile);
      
      // Assert
      expect(result).toBeSuccessResponse();
      expect(result.data.chemicals).toContainMolecules(1);
      expect(result.data.chemicals[0]).toBeMolecularData();
    });
    
    it('should reject oversized images', async () => {
      // Arrange
      const largeImage = create.file()
        .asPng()
        .asLargeFile() // 10MB
        .build();
      
      // Act & Assert
      await expect(processor.processImage(largeImage))
        .rejects.toThrow('Image size exceeds maximum allowed');
    });
    
    it('should validate image format', async () => {
      // Arrange
      const invalidFile = create.file()
        .withMimetype('text/plain')
        .withName('not-an-image.txt')
        .build();
      
      // Act & Assert
      await expect(processor.processImage(invalidFile))
        .rejects.toThrow('Invalid image format');
    });
  });
  
  describe('generateSDF', () => {
    it('should generate SDF file for molecule', async () => {
      // Arrange
      const molecule = create.molecule()
        .withName('Ethanol')
        .withSmiles('CCO')
        .build();
      
      const fs = require('fs');
      fs.writeFileSync = jest.fn();
      fs.existsSync = jest.fn().mockReturnValue(false);
      fs.mkdirSync = jest.fn();
      
      // Act
      const result = await processor.generateSDF(molecule);
      
      // Assert
      expect(result).toMatchShape({
        success: 'boolean',
        path: 'string',
        molecule: 'object'
      });
      expect(result.path).toMatch(/ethanol.*\.sdf$/);
      expect(fs.writeFileSync).toHaveBeenCalled();
    });
    
    it('should sanitize molecule names for filenames', async () => {
      // Arrange
      const molecule = create.molecule()
        .withName('Molecule/With\\Special*Chars')
        .withSmiles('CCO')
        .build();
      
      const fs = require('fs');
      fs.writeFileSync = jest.fn();
      fs.existsSync = jest.fn().mockReturnValue(false);
      fs.mkdirSync = jest.fn();
      
      // Act
      const result = await processor.generateSDF(molecule);
      
      // Assert
      expect(result.path).not.toMatch(/[\/\\*]/);
      expect(result.path).toMatch(/molecule-with-special-chars.*\.sdf$/);
    });
  });
  
  describe('batchProcess', () => {
    it('should process multiple items in parallel', async () => {
      // Arrange
      const items = [
        { type: 'text', content: 'What is in water?' },
        { type: 'text', content: 'What is in ethanol?' },
        { type: 'text', content: 'What is in caffeine?' }
      ];
      
      mockOpenAI.mockImplementation(() => ({
        chat: {
          completions: {
            create: jest.fn()
              .mockResolvedValueOnce({
                choices: [{
                  message: {
                    content: JSON.stringify({
                      object: 'water',
                      chemicals: [create.molecule().asWater().build()]
                    })
                  }
                }]
              })
              .mockResolvedValueOnce({
                choices: [{
                  message: {
                    content: JSON.stringify({
                      object: 'ethanol',
                      chemicals: [create.molecule().build()]
                    })
                  }
                }]
              })
              .mockResolvedValueOnce({
                choices: [{
                  message: {
                    content: JSON.stringify({
                      object: 'caffeine',
                      chemicals: [create.molecule().asCaffeine().build()]
                    })
                  }
                }]
              })
          }
        }
      }));
      
      // Act
      const results = await processor.batchProcess(items);
      
      // Assert
      expect(results).toHaveLength(3);
      results.forEach(result => {
        expect(result).toBeSuccessResponse();
        expect(result.data.chemicals).toContainMolecules();
      });
      
      // Verify parallel execution
      const mockCreate = mockOpenAI.mock.results[0].value.chat.completions.create;
      expect(mockCreate).toHaveBeenCalledTimes(3);
    });
    
    it('should handle partial failures in batch', async () => {
      // Arrange
      const items = [
        { type: 'text', content: 'Valid request' },
        { type: 'text', content: '' }, // Invalid
        { type: 'text', content: 'Another valid request' }
      ];
      
      mockOpenAI.mockImplementation(() => ({
        chat: {
          completions: {
            create: jest.fn()
              .mockResolvedValueOnce({
                choices: [{
                  message: {
                    content: JSON.stringify({
                      object: 'test',
                      chemicals: [create.molecule().build()]
                    })
                  }
                }]
              })
              .mockResolvedValueOnce({
                choices: [{
                  message: {
                    content: JSON.stringify({
                      object: 'test2',
                      chemicals: [create.molecule().build()]
                    })
                  }
                }]
              })
          }
        }
      }));
      
      // Act
      const results = await processor.batchProcess(items);
      
      // Assert
      expect(results).toHaveLength(3);
      expect(results[0]).toBeSuccessResponse();
      expect(results[1]).toBeErrorResponse();
      expect(results[2]).toBeSuccessResponse();
    });
  });
});
