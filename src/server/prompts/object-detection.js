function buildObjectDetectionPrompt(x, y) {
  const focus = (typeof x === 'number' && typeof y === 'number')
    ? `Focus near (${Math.round(x)}, ${Math.round(y)}).`
    : '';
  return [
    'Task: Identify the real-world object in the image and return JSON only.',
    focus,
    'Return the most specific object description (e.g., "red apple", "PET plastic bottle").',
    'If possible, suggest a bounding box around the object of interest.',
    'JSON format:',
    '{ "object": "string", "recommendedBox": { "x": number, "y": number, "width": number, "height": number } }'
  ].join('\n');
}

module.exports = { buildObjectDetectionPrompt };


