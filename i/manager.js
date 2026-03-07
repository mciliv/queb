#!/usr/bin/env node

require('dotenv').config();
const { spawn } = require('child_process');
const path = require('path');
const { aiManager } = require('./AIManager.cjs');
const { report, fetchAgents, sendDirective, getId } = require('./lib/apm/client');
const { start: startServer } = require('./lib/apm/server');

const children = new Map();

const SYSTEM = `You are a coding manager. You oversee multiple AI coding agents working in parallel.

You can see what each agent is doing — their status, current task, which files they're touching, errors they hit.
Since they're language models like you, you can read their prompts and understand their intent semantically.

Your job:
- Break down the user's goal into parallel tasks when possible
- Spawn agents for each task
- Watch their progress and intervene if they get stuck or conflict
- Send directives to redirect agents
- Report back when the goal is complete

Keep spawned tasks focused and small. Don't spawn more than 4 agents at once.
When all tasks complete, summarize what was done.`;

const tools = {
  check_agents: {
    description: 'See all active agents, their status, what they are working on, errors, etc.',
    parameters: {},
    execute: async () => {
      const agents = await fetchAgents();
      if (!agents.length) return 'No agents currently active.';
      return agents.map(a => {
        const age = Math.round((Date.now() - a.lastUpdate) / 1000);
        let line = `[${a.id.slice(0, 8)}] status=${a.status} name=${a.name || '?'}`;
        if (a.agent) line += ` agent=${a.agent}`;
        if (a.provider) line += ` provider=${a.provider}/${a.model || '?'}`;
        if (a.prompt) line += `\n  prompt: "${a.prompt}"`;
        if (a.tool) line += `\n  tool: ${a.tool}`;
        if (a.error) line += `\n  ERROR: ${a.error}`;
        if (a.durationMs) line += `\n  duration: ${a.durationMs}ms`;
        line += `\n  last seen: ${age}s ago`;
        return line;
      }).join('\n\n');
    }
  },

  spawn_agent: {
    description: 'Spawn a new coding agent with a specific task. The agent runs cli.js which can read/write files and run commands.',
    parameters: { task: 'string - the prompt/task for the agent' },
    execute: async ({ task }) => {
      const cliPath = path.join(__dirname, 'cli.js');
      const child = spawn('node', [cliPath, task], {
        cwd: process.cwd(),
        stdio: ['ignore', 'pipe', 'pipe'],
        env: { ...process.env }
      });

      let output = '';
      child.stdout.on('data', d => output += d);
      child.stderr.on('data', d => output += d);

      // Wait briefly for it to register with APM
      await new Promise(r => setTimeout(r, 500));

      // Find its agent ID from APM
      const agents = await fetchAgents();
      const recent = agents
        .filter(a => a.status === 'running' && a.prompt?.includes(task.slice(0, 30)))
        .sort((a, b) => b.lastUpdate - a.lastUpdate)[0];

      const agentId = recent?.id || 'unknown';
      children.set(agentId, { child, task, output: () => output });

      child.on('exit', (code) => {
        report({ name: 'manager', status: 'running', lastEvent: `agent ${agentId.slice(0, 8)} exited (code ${code})` });
      });

      return `Spawned agent ${agentId.slice(0, 8)} for: "${task.slice(0, 100)}"`;
    }
  },

  message_agent: {
    description: 'Send a directive/message to a specific agent (by ID prefix). The agent will see this on its next check.',
    parameters: { agent_id: 'string - agent ID or prefix', message: 'string - the directive' },
    execute: async ({ agent_id, message }) => {
      const agents = await fetchAgents();
      const match = agents.find(a => a.id.startsWith(agent_id));
      if (!match) return `No agent found matching "${agent_id}"`;
      sendDirective(match.id, message);
      return `Directive sent to ${match.id.slice(0, 8)}: "${message.slice(0, 100)}"`;
    }
  },

  read_agent_output: {
    description: 'Read the stdout/stderr output of a spawned agent.',
    parameters: { agent_id: 'string - agent ID or prefix' },
    execute: async ({ agent_id }) => {
      for (const [id, info] of children) {
        if (id.startsWith(agent_id)) {
          const out = info.output();
          return out.slice(-3000) || '(no output yet)';
        }
      }
      return `No spawned agent found matching "${agent_id}"`;
    }
  },

  wait: {
    description: 'Wait for a number of seconds, then check agent status. Use when agents are working and you need to give them time.',
    parameters: { seconds: 'number - how long to wait (max 30)' },
    execute: async ({ seconds }) => {
      const secs = Math.min(Number(seconds) || 5, 30);
      await new Promise(r => setTimeout(r, secs * 1000));
      // Return current status after waiting
      const agents = await fetchAgents();
      if (!agents.length) return `Waited ${secs}s. No agents active.`;
      return `Waited ${secs}s. ${agents.length} agent(s): ` +
        agents.map(a => `${a.name || a.id.slice(0, 8)}=${a.status}`).join(', ');
    }
  }
};

async function run(goal) {
  report({ name: 'manager', status: 'running', agent: 'manager', prompt: goal.slice(0, 200) });
  console.log(`\x1b[36m[Manager]\x1b[0m Goal: ${goal}\n`);

  const messages = [
    { role: 'system', content: SYSTEM },
    { role: 'user', content: goal }
  ];

  for (let i = 0; i < 30; i++) {
    try {
      const response = await aiManager.execute('planner', { messages, tools });

      if (response.toolCalls?.length) {
        messages.push({ role: 'assistant', content: response.text || '', toolCalls: response.toolCalls });

        for (const tc of response.toolCalls) {
          console.log(`\x1b[35m[${tc.toolName}]\x1b[0m`, tc.args?.task || tc.args?.message || tc.args?.agent_id || '');

          const tool = tools[tc.toolName];
          const result = tool ? await tool.execute(tc.args || {}) : `Unknown tool: ${tc.toolName}`;
          console.log(`\x1b[90m${result.slice(0, 200)}\x1b[0m\n`);

          messages.push({ role: 'tool', content: [{ type: 'tool-result', toolName: tc.toolName, toolCallId: tc.toolCallId, result }] });
        }
      } else {
        console.log(`\n\x1b[32m[Manager Done]\x1b[0m\n${response.text}`);
        report({ name: 'manager', status: 'idle' });
        break;
      }
    } catch (error) {
      console.error(`\x1b[31m[Manager Error]\x1b[0m`, error.message);
      report({ name: 'manager', status: 'error', error: error.message });
      break;
    }
  }

  // Cleanup children
  for (const [, info] of children) {
    if (!info.child.killed) info.child.kill();
  }
}

// Entry point
if (require.main === module) {
  const args = process.argv.slice(2);

  if (args[0] === '--server' || args[0] === '-s') {
    // Start APM server + manager together
    startServer();
    const goal = args.slice(1).join(' ');
    if (goal) setTimeout(() => run(goal), 500);
    else console.log('APM server started. Pass a goal to start managing agents.');
  } else {
    const goal = args.join(' ');
    if (!goal) {
      console.log('Usage: node i/manager.js [-s] "goal"');
      console.log('  -s  Also start the APM server');
      console.log('\nExample: node i/manager.js -s "add input validation to api routes and write tests"');
      process.exit(1);
    }
    run(goal);
  }
}

module.exports = { run };
