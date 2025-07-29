/**
 * ATOM PREDICTOR SERVICE
 * Purpose: OpenAI-powered molecular analysis and chemical identification
 * Features: Image analysis, text analysis, SMILES generation, structured chemical data
 * 
 * Core Responsibilities:
 * - Analyze images to identify chemical compounds
 * - Process text descriptions to extract molecular information
 * - Generate accurate SMILES notation for identified chemicals
 * - Provide structured chemical data for 3D visualization
 */

const { OpenAI } = require("openai");
const {
  ObjectIdentificationSchema,
  CHEMICAL_REPRESENTATIONS,
} = require("../schemas/schemas");
const ErrorHandler = require("./error-handler");
const { AI_CONFIG } = require("../config/constants");

class AtomPredictor {
  constructor(apiKey) {
    this.client = new OpenAI({ 
      apiKey,
      timeout: AI_CONFIG.AI_TIMEOUT,
      maxRetries: AI_CONFIG.MAX_RETRIES
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

CRITICAL SMILES Guidelines:
- Use ONLY standard SMILES notation, never molecular formulas
- Provide accurate SMILES for the actual molecules present
- Examples: "O" (water), "CCO" (ethanol), "C1=CC=CC=C1" (benzene)
- For complex molecules, provide the complete accurate SMILES

NEVER use molecular formulas like "H2O", "C2H6O", "CaCO3" - always use SMILES.
Ensure SMILES are chemically accurate and can be parsed by standard tools.`;
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
            detail: AI_CONFIG.IMAGE_DETAIL,
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
        model: AI_CONFIG.MODEL,
        messages,
        max_tokens: AI_CONFIG.MAX_TOKENS,
        temperature: AI_CONFIG.TEMPERATURE,
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
        model: AI_CONFIG.MODEL,
        messages: [
          {
            role: "user",
            content: `Analyze this object: "${object}". ${this.chemicalInstructions}`,
          },
        ],
        max_tokens: AI_CONFIG.MAX_TOKENS,
        temperature: AI_CONFIG.TEMPERATURE,
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
