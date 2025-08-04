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
1. BIOLOGICAL MATERIALS (skin, food, plants, etc.): Provide major molecular components typically found in that material
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

      
      const startTime = Date.now();
      
      const response = await this.client.chat.completions.create({
        model: "gpt-4o",
        messages,
        max_tokens: 1000,
        temperature: 0.1,
      });

      const requestTime = Date.now() - startTime;
      

      const content = response.choices[0].message.content;
      
      const parsed = this.parseAIResponse(content);
      

      return {
        object: parsed.object || "Unknown object",
        chemicals: parsed.chemicals || [],
      };
    } catch (error) {
      console.error('💥 AtomPredictor error:', error);
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
    try {
      // Try to extract JSON from the response
      let jsonMatch = content.match(/\{[\s\S]*\}/);
      if (jsonMatch) {
        let jsonStr = jsonMatch[0];
        
        // If JSON appears truncated, try to fix it
        if (!jsonStr.endsWith('}')) {
          // Find the last complete chemical entry
          const chemicalsMatch = jsonStr.match(/"chemicals"\s*:\s*\[([\s\S]*)/);
          if (chemicalsMatch) {
            const chemicalsStr = chemicalsMatch[1];
            const lastCompleteEntry = chemicalsStr.lastIndexOf('{"name"');
            if (lastCompleteEntry > 0) {
              const truncatedPart = chemicalsStr.substring(0, lastCompleteEntry - 1);
              jsonStr = jsonStr.replace(/"chemicals"\s*:\s*\[[\s\S]*/, `"chemicals": [${truncatedPart}]}`);
            } else {
              // No valid chemicals found, return empty array
              const objectMatch = jsonStr.match(/"object"\s*:\s*"([^"]+)"/);
              const objectName = objectMatch ? objectMatch[1] : "Unknown object";
              jsonStr = `{"object": "${objectName}", "chemicals": []}`;
            }
          }
        }
        
        return JSON.parse(jsonStr);
      }

      // Fallback: try to parse the entire content as JSON
      return JSON.parse(content);
    } catch (error) {
      console.error("Failed to parse AI response:", content.substring(0, 500) + "...");
      
      // Try to extract object name at least
      const objectMatch = content.match(/"object"\s*:\s*"([^"]+)"/);
      const objectName = objectMatch ? objectMatch[1] : "Unknown object";
      
      return {
        object: objectName,
        chemicals: []
      };
    }
  }
}

module.exports = AtomPredictor;
