const { report } = require('./apm/client');

class AgentRegistry {
  constructor(modelManager) {
    this.modelManager = modelManager;
    this.agents = {
      coder: {
        role: 'Senior Software Engineer',
        priorityList: ['openai', 'xai', 'google'],
        instructions: 'Expert at Node.js, Python, and system architecture. Focus on clean, efficient, and testable code.'
      },
      architect: {
        role: 'Systems Architect',
        priorityList: ['google', 'openai', 'xai'],
        instructions: 'Expert at designing scalable systems. Focus on modularity, security, and performance.'
      },
      researcher: {
        role: 'Technical Researcher',
        priorityList: ['xai', 'openai', 'google'],
        instructions: 'Expert at finding and synthesizing technical information.'
      }
    };
  }

  async runTask(agentName, prompt, params = {}) {
    const agent = this.agents[agentName];
    if (!agent) throw new Error(`Unknown agent: ${agentName}`);

    const start = Date.now();
    report({ name: agentName, status: 'running', agent: agentName, prompt: prompt?.slice(0, 200) });

    const taskParams = {
      ...params,
      messages: [
        { role: 'system', content: `You are a ${agent.role}. ${agent.instructions}` },
        ...(params.messages || []),
        { role: 'user', content: prompt }
      ]
    };

    try {
      const result = await this.modelManager.callWithFailover(taskParams, agent.priorityList);
      report({ name: agentName, status: 'idle', agent: agentName, durationMs: Date.now() - start });
      return result;
    } catch (error) {
      report({ name: agentName, status: 'error', agent: agentName, error: error.message });
      throw error;
    }
  }

  getAgents() {
    return Object.keys(this.agents);
  }
}

module.exports = AgentRegistry;
