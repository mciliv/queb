/**
 * Queb Client Hooks
 * 
 * These hooks provide clear, user-focused interfaces for different
 * aspects of molecular analysis. Each hook is designed around a
 * specific user action or goal.
 */

// Analysis hooks - "I want to analyze..."
export { useAnalyzeText } from './useAnalyzeText';       // ...text descriptions
export { useAnalyzeImage } from './useAnalyzeImage';     // ...uploaded images  
export { useAnalyzeCamera } from './useAnalyzeCamera';   // ...live camera feed

// Data hooks - "I want to know about..."
export { useChemicalData } from './useChemicalData';     // ...a specific chemical
export { useMolecule3D } from './useMolecule3D';         // ...3D structures

// Original hook (still available for backward compatibility)
export { useApi } from './useApi';

/**
 * Usage Examples:
 * 
 * // Analyzing text
 * const { analyze, molecules } = useAnalyzeText();
 * await analyze('coffee');
 * 
 * // Analyzing images
 * const { analyzeImage, molecules } = useAnalyzeImage();
 * await analyzeImage(imageFile);
 * 
 * // Camera analysis
 * const { startCamera, analyzeClick } = useAnalyzeCamera();
 * await startCamera(videoElement);
 * await analyzeClick(x, y);
 * 
 * // Looking up chemical data
 * const { lookupChemical, chemicalInfo } = useChemicalData();
 * await lookupChemical('caffeine');
 * 
 * // Generating 3D structures
 * const { generate3DFromName } = useMolecule3D();
 * await generate3DFromName('aspirin');
 */


