#!/usr/bin/env node

require('dotenv').config();
const fs = require('fs');
const path = require('path');
const { execSync } = require('child_process');
const readline = require('readline');
const { z } = require('zod');
const { AIService } = require('./index');
const { report, getDirectives } = require('./lib/apm/client');

// --- Configuration ---
const MAX_ITERATIONS = 20;
const DEBUG = process.env.DEBUG === 'true';

// --- Tools Definition ---
const tools = {
  read_file: {
    description: 'Read the content of a file',
    parameters: z.object({
      file_path: z.string().describe('The path to the file to read'),
    }),
    execute: async ({ file_path }) => {
      try {
        const content = fs.readFileSync(file_path, 'utf8');
        return `[read_file] Content of ${file_path}:
${content.slice(0, 10000)}${content.length > 10000 ? '
...(truncated)...' : ''}`;
      } catch (error) {
        return `[read_file] Error: ${error.message}`;
      }
    },
  },
  write_file: {
    description: 'Write content to a file',
    parameters: z.object({
      file_path: z.string().describe('The path to the file to write'),
      content: z.string().describe('The content to write'),
    }),
    execute: async ({ file_path, content }) => {
      try {
        const dir = path.dirname(file_path);
        if (!fs.existsSync(dir)) fs.mkdirSync(dir, { recursive: true });
        fs.writeFileSync(file_path, content, 'utf8');
        return `[write_file] Successfully wrote to ${file_path}`;
      } catch (error) {
        return `[write_file] Error: ${error.message}`;
      }
    },
  },
  run_shell_command: {
    description: 'Execute a shell command',
    parameters: z.object({
      command: z.string().describe('The command to execute'),
    }),
    execute: async ({ command }) => {
      console.log(`\x1b[33m[EXEC]\x1b[0m ${command}`);
      // Safety check could go here
      try {
        const output = execSync(command, { encoding: 'utf8', stdio: 'pipe' });
        return `[run_shell_command] Output:
${output}`;
      } catch (error) {
        return `[run_shell_command] Error (${error.code}):
${error.stdout}
${error.stderr}`;
      }
    },
  },
  list_files: {
    description: 'List files in a directory',
    parameters: z.object({
      dir_path: z.string().describe('Directory path to list'),
      recursive: z.boolean().optional().describe('Recursive listing'),
    }),
    execute: async ({ dir_path, recursive }) => {
      try {
        // Simple implementation using find or ls
        const cmd = recursive ? `find "${dir_path}" -maxdepth 3 -not -path '*/.*'` : `ls -F "${dir_path}"`;
        const output = execSync(cmd, { encoding: 'utf8' });
        return `[list_files] Result:
${output}`;
      } catch (error) {
        return `[list_files] Error: ${error.message}`;
      }
    },
  }
};

// --- Agent Logic ---

async function runAgent(prompt) {
  const ai = new AIService();
  const start = Date.now();

  console.log(`\x1b[36m[Agent]\x1b[0m Started with provider: ${ai.getProvider()} model: ${ai.getModel()}`);
  report({ name: 'cli', status: 'running', agent: 'cli', prompt: prompt?.slice(0, 200), provider: ai.getProvider(), model: ai.getModel() });

  let messages = [
    { role: 'system', content: 'You are a helpful AI coding assistant. You have access to tools to read files, write files, and run commands. Use them to fulfill the user request. When you are done, provide a final answer.' },
    { role: 'user', content: prompt }
  ];

  let iterations = 0;

  while (iterations < MAX_ITERATIONS) {
    iterations++;
    if (DEBUG) console.log(`\n--- Iteration ${iterations} ---`);

    try {
      // Check for directives from manager
      const directives = await getDirectives();
      for (const d of directives) {
        console.log(`\x1b[33m[Directive]\x1b[0m ${d.message}`);
        messages.push({ role: 'user', content: `[Manager directive]: ${d.message}` });
      }

      const response = await ai.callAPI({ messages, tools });

      if (response.toolCalls && response.toolCalls.length > 0) {
        messages.push({ role: 'assistant', content: response.content || '', toolCalls: response.toolCalls });

        for (const toolCall of response.toolCalls) {
          console.log(`\x1b[35m[Tool]\x1b[0m Calling ${toolCall.toolName}...`);
          report({ name: 'cli', status: 'running', tool: toolCall.toolName });

          const result = tools[toolCall.toolName]
            ? await tools[toolCall.toolName].execute(toolCall.args)
            : `Error: Tool ${toolCall.toolName} not found`;

          messages.push({ role: 'tool', content: [{ type: 'tool-result', toolName: toolCall.toolName, toolCallId: toolCall.toolCallId, result }] });
        }
      } else {
        console.log(`\n\x1b[32m[Result]\x1b[0m\n${response.content}`);
        report({ name: 'cli', status: 'idle', durationMs: Date.now() - start });
        return;
      }
    } catch (error) {
      console.error(`\x1b[31m[Error]\x1b[0m`, error);
      report({ name: 'cli', status: 'error', error: error.message });
      break;
    }
  }
}

// --- Entry Point ---
if (require.main === module) {
  const args = process.argv.slice(2);
  const prompt = args.join(' ');

  if (!prompt) {
    console.log('Usage: node cli.js <prompt>');
    process.exit(1);
  }

  runAgent(prompt);
}

module.exports = { runAgent };
