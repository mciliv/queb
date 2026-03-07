#!/usr/bin/env node

import { xai } from '@ai-sdk/xai';
import { generateText } from 'ai';

export async function query(system, prompt) {
  const { text, response } = await generateText({
    model: xai.responses('grok-4'),
    system: system,
    prompt: prompt,
  });

  console.log(text);

  // The response ID can be used to continue the conversation
  console.log(response.id);
  
  return { text, response };
}

if (import.meta.url === `file://${process.argv[1]}`) {
  const [system, prompt] = process.argv.slice(2);
  if (!system || !prompt) {
    console.error('Usage: node x.js "<system>" "<prompt>"');
    process.exit(1);
  }
  query(system, prompt).catch((error) => {
    console.error(error?.message || String(error));
    process.exit(1);
  });
}
