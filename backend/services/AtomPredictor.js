const { OpenAI } = require("openai");
const {
  ObjectIdentificationSchema,
  CHEMICAL_REPRESENTATIONS,
} = require("../schemas/schemas");

// Import the new prompt engineering modules
const { buildChemicalAnalysisInstructions } = require("../prompts/chemical-analysis-instructions");
const { parseAIResponseWithFallbacks, validateSMILESQuality } = require("../prompts/fallback-handlers");
const { getRelevantExamples } = require("../prompts/material-examples");

class AtomPredictor {
  constructor(apiKey) {
    this.client = new OpenAI({ apiKey });
    // Use the improved instruction builder from git history analysis
    this.chemicalInstructions = this.buildChemicalInstructions();
  }

  buildChemicalInstructions() {
    // Use the comprehensive instructions with proven techniques
    return buildChemicalAnalysisInstructions();
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

      // Add cropped region if available - improved focus text from git history
      if (croppedImageBase64) {
        let focusText = `Here's a cropped view of the area of interest. Analyze the chemical composition of the material or substance visible in this region. Use the examples above as your guide for accurate SMILES notation:`;

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
        temperature: 0.1, // Keep low for consistency
      });

      const content = response.choices[0].message.content;
      
      // Use the improved parsing with smart fallbacks
      const parsed = parseAIResponseWithFallbacks(content);
      
      // Validate and improve SMILES quality
      const validatedChemicals = validateSMILESQuality(parsed.chemicals || []);

      return {
        object: parsed.object || "Unknown object",
        chemicals: validatedChemicals,
      };
    } catch (error) {
      console.error("AI analysis error:", error);
      throw new Error(`AI analysis failed: ${error.message}`);
    }
  }

  async analyzeText(object) {
    try {
      // Enhanced text analysis with context-aware examples
      const objectType = this.detectObjectType(object);
      const relevantExamples = getRelevantExamples(objectType);
      
      const enhancedInstructions = `${this.chemicalInstructions}

${relevantExamples}

Now analyze this specific object: "${object}"`;

      const response = await this.client.chat.completions.create({
        model: "gpt-4o",
        messages: [
          {
            role: "user",
            content: enhancedInstructions,
          },
        ],
        max_tokens: 1000,
        temperature: 0.1,
      });

      const content = response.choices[0].message.content;
      
      // Use the improved parsing with smart fallbacks
      const parsed = parseAIResponseWithFallbacks(content);
      
      // Validate and improve SMILES quality
      const validatedChemicals = validateSMILESQuality(parsed.chemicals || []);

      return {
        object: parsed.object || object,
        chemicals: validatedChemicals,
      };
    } catch (error) {
      console.error("AI text analysis error:", error);
      throw new Error(`AI text analysis failed: ${error.message}`);
    }
  }

  /**
   * Detect object type for context-aware examples
   * Helps provide more relevant prompts to improve accuracy
   */
  detectObjectType(object) {
    const objectLower = object.toLowerCase();
    
    if (objectLower.includes('wine') || objectLower.includes('beer') || 
        objectLower.includes('coffee') || objectLower.includes('drink') ||
        objectLower.includes('beverage')) {
      return 'beverage';
    }
    
    if (objectLower.includes('fruit') || objectLower.includes('food') ||
        objectLower.includes('apple') || objectLower.includes('orange') ||
        objectLower.includes('vegetable')) {
      return 'food';
    }
    
    if (objectLower.includes('plastic') || objectLower.includes('metal') ||
        objectLower.includes('wood') || objectLower.includes('stone')) {
      return 'material';
    }
    
    return 'general';
  }

  // Keep the legacy parseAIResponse method for compatibility but mark as deprecated
  parseAIResponse(content) {
    console.warn("Using deprecated parseAIResponse - consider updating to use fallback-handlers module");
    return parseAIResponseWithFallbacks(content);
  }
}

module.exports = AtomPredictor;
