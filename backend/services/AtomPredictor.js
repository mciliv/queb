const { OpenAI } = require("openai");
const {
  ObjectIdentificationSchema,
  CHEMICAL_REPRESENTATIONS,
} = require("../schemas/schemas");
const ErrorHandler = require("./error-handler");

class AtomPredictor {
  constructor(apiKey) {
    this.client = new OpenAI({ 
      apiKey,
      timeout: 30000, // 30 seconds timeout
      maxRetries: 2   // Retry failed requests twice
    });
    this.chemicalInstructions = this.buildChemicalInstructions();
  }

  buildChemicalInstructions() {
    return `Analyze the object and provide a JSON response with chemical components.

Response format:
{
  "object": "Object name",
  "chemicals": [
    {"name": "Chemical name", "smiles": "SMILES notation"}
  ]
}

ANALYSIS APPROACH:
1. BIOLOGICAL MATERIALS (skin, food, plants, etc.): Provide all molecular components typically found in that material
2. MANUFACTURED ITEMS: Analyze visible materials and compositions
3. PURE SUBSTANCES: Identify the specific chemical if recognizable

BIOLOGICAL MATERIAL EXAMPLES:
- Human skin: Include keratin proteins, lipids, water, cholesterol, ceramides
- Food items: Include major nutrients, flavor compounds, preservatives
- Plants: Include cellulose, chlorophyll, common plant metabolites

CRITICAL SMILES Guidelines:
- Use ONLY standard SMILES notation, never molecular formulas
- Provide accurate SMILES for the actual molecules present
- Examples: "O" (water), "CCO" (ethanol), "C1=CC=CC=C1" (benzene)
- For proteins, use representative amino acids: "N[C@@H](CC1=CC=CC=C1)C(=O)O" (phenylalanine)
- For complex molecules, provide the complete accurate SMILES

NEVER use molecular formulas like "H2O", "C2H6O", "CaCO3" - always use SMILES.
When visual analysis is insufficient, use scientific knowledge of typical molecular composition.`;
  }

  async analyzeImage(
    imageBase64,
    croppedImageBase64 = null,
    x = null,
    y = null,
    cropMiddleX = null,
    cropMiddleY = null,
    cropSize = null,
  ) {
    try {
      const messages = [
        {
          role: "user",
          content: [
            {
              type: "text",
              text: this.chemicalInstructions,
            },
            {
              type: "image_url",
              image_url: {
                url: `data:image/jpeg;base64,${imageBase64}`,
                detail: "high",
              },
            },
          ],
        },
      ];

      // Add cropped region if available
      if (croppedImageBase64) {
        let focusText = `Here's a cropped view of the area of interest. Analyze the chemical composition of the material or substance visible in this region:`;

        messages[0].content.push({
          type: "text",
          text: focusText,
        });
        messages[0].content.push({
          type: "image_url",
          image_url: {
            url: `data:image/jpeg;base64,${croppedImageBase64}`,
            detail: "high",
          },
        });
      }

      const response = await this.client.chat.completions.create({
        model: "gpt-4o",
        messages,
        max_tokens: 1000,
        temperature: 0.1,
      });

      const content = response.choices[0].message.content;
      const parsed = this.parseAIResponse(content);

      return {
        object: parsed.object || "Unknown object",
        chemicals: parsed.chemicals || [],
      };
    } catch (error) {
      const { errorMessage } = ErrorHandler.handleAIError(error, 'image analysis');
      throw new Error(errorMessage);
    }
  }

  async analyzeText(object) {
    try {
      const response = await this.client.chat.completions.create({
        model: "gpt-4o",
        messages: [
          {
            role: "user",
            content: `Analyze this object: "${object}". ${this.chemicalInstructions}`,
          },
        ],
        max_tokens: 1000,
        temperature: 0.1,
      });

      const content = response.choices[0].message.content;
      const parsed = this.parseAIResponse(content);

      return {
        object: parsed.object || object,
        chemicals: parsed.chemicals || [],
      };
    } catch (error) {
      const { errorMessage } = ErrorHandler.handleAIError(error, 'text analysis');
      throw new Error(errorMessage);
    }
  }

  parseAIResponse(content) {
    // Handle null/undefined content first
    if (!content) {
      return {
        object: "Unknown object",
        chemicals: []
      };
    }

    try {
      // Extract JSON from potential markdown or text wrapper
      const jsonMatch = content.match(/```(?:json)?\s*(\{[\s\S]*?\})\s*```/);
      if (jsonMatch) {
        return JSON.parse(jsonMatch[1]);
      }

      // Look for standalone JSON object in the response
      const objectMatch = content.match(/\{[\s\S]*?\}/);
      if (objectMatch) {
        const jsonStr = objectMatch[0];
        
        // Check if this appears to be a complete, properly structured response
        if (jsonStr.includes('"object"') && jsonStr.includes('"chemicals"')) {
          // Validate that all chemical entries are complete
          const chemicalsMatch = jsonStr.match(/"chemicals"\s*:\s*\[([\s\S]*?)\]/);
          if (chemicalsMatch) {
            const chemicalsStr = chemicalsMatch[1];
            const entryCount = (chemicalsStr.match(/\{/g) || []).length;
            const completeEntries = (chemicalsStr.match(/\{\s*"name"\s*:\s*"[^"]*"\s*,\s*"smiles"\s*:\s*"[^"]*"\s*\}/g) || []).length;
            
            // If some entries are incomplete, extract what we can
            if (entryCount > completeEntries) {
              const completeChemicals = [];
              const completeMatches = chemicalsStr.match(/\{\s*"name"\s*:\s*"([^"]*)"\s*,\s*"smiles"\s*:\s*"([^"]*)"\s*\}/g) || [];
              
              completeMatches.forEach(match => {
                const nameMatch = match.match(/"name"\s*:\s*"([^"]*)"/);
                const smilesMatch = match.match(/"smiles"\s*:\s*"([^"]*)"/);
                if (nameMatch && smilesMatch) {
                  completeChemicals.push({
                    name: nameMatch[1],
                    smiles: smilesMatch[1]
                  });
                }
              });
              
              // Get object name
              const objectMatch = jsonStr.match(/"object"\s*:\s*"([^"]+)"/);
              const objectName = objectMatch ? objectMatch[1] : "Unknown object";
              
              return {
                object: objectName,
                chemicals: completeChemicals
              };
            }
          }
        }
        
        return JSON.parse(jsonStr);
      }

      // Fallback: try to parse the entire content as JSON
      return JSON.parse(content);
    } catch (error) {
      // Enhanced fallback: try to extract chemicals even from malformed JSON
      const objectMatch = content.match(/"object"\s*:\s*"([^"]+)"/);
      const objectName = objectMatch ? objectMatch[1] : "Unknown object";
      
      // Try to extract individual chemical entries
      const chemicals = [];
      const entryRegex = /\{"name"\s*:\s*"([^"]+)"\s*,\s*"smiles"\s*:\s*"([^"]+)"\}/g;
      let match;
      
      while ((match = entryRegex.exec(content)) !== null) {
        chemicals.push({
          name: match[1],
          smiles: match[2]
        });
      }
      
      return {
        object: objectName,
        chemicals: chemicals
      };
    }
  }
}

module.exports = AtomPredictor;
