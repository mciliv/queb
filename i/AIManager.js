/**
 * AIManager: Standalone AI Orchestrator
 * - Handles direct SDK calls (OpenAI, xAI, Google)
 * - Implements Model Priority Queue for Rate Limit (429) Failover
 * - Supports Specialized Agents for different tasks
 */

import { openai } from '@ai-sdk/openai';
import { xai } from '@ai-sdk/xai';
import { google } from '@ai-sdk/google';
import { generateText } from 'ai';

class AIManager {
  constructor() {
    this.throttledProviders = new Map();
    this.backoffTime = 60000; // 1 minute

    // Define Specialized Agents & their Priority Queues
    this.agents = {
      coder: {
        description: 'Optimized for high-fidelity code generation',
        priority: ['openai', 'xai', 'google'],
        models: {
          openai: 'gpt-4o',
          xai: 'grok-beta',
          google: 'gemini-1.5-pro'
        }
      },
      planner: {
        description: 'Optimized for massive context and reasoning',
        priority: ['google', 'openai', 'xai'],
        models: {
          google: 'gemini-2.0-flash',
          openai: 'o1-preview',
          xai: 'grok-beta'
        }
      },
      fast: {
        description: 'Low-latency small model fallback',
        priority: ['google', 'openai'],
        models: {
          google: 'gemini-1.5-flash',
          openai: 'gpt-4o-mini'
        }
      }
    };
  }

  /**
   * Get the SDK model object for a given provider
   */
  _getModel(providerName, agentName) {
    const modelId = this.agents[agentName].models[providerName];
    
    switch (providerName) {
      case 'openai': return openai(modelId);
      case 'xai': return xai(modelId);
      case 'google': return google(modelId);
      default: throw new Error(`Unsupported provider: ${providerName}`);
    }
  }

  /**
   * Execute task with automatic failover
   */
  async execute(agentName, params) {
    const agent = this.agents[agentName];
    if (!agent) throw new Error(`Unknown agent: ${agentName}`);

    const priorityList = agent.priority;
    let lastError = null;

    for (const provider of priorityList) {
      if (this._isThrottled(provider)) {
        console.warn(`[AIManager] Skipping throttled provider: ${provider}`);
        continue;
      }

      try {
        const model = this._getModel(provider, agentName);
        console.log(`[AIManager] Attempting task with ${provider} (${agent.models[provider]})`);

        const result = await generateText({
          model: model,
          ...params
        });

        return {
          ...result,
          providerUsed: provider,
          modelUsed: agent.models[provider]
        };

      } catch (error) {
        lastError = error;
        
        if (this._isRateLimitError(error)) {
          console.error(`[AIManager] 429 Rate Limit hit for ${provider}. Failing over...`);
          this._markThrottled(provider);
          continue; 
        }

        console.error(`[AIManager] Provider ${provider} failed: ${error.message}`);
        // For non-rate-limit errors, we still attempt failover unless it's a critical logic error
        continue;
      }
    }

    throw new Error(`AIManager: All providers in queue failed. Last error: ${lastError?.message}`);
  }

  _isRateLimitError(error) {
    const msg = error.message.toLowerCase();
    return msg.includes('429') || msg.includes('rate limit') || msg.includes('too many requests');
  }

  _markThrottled(provider) {
    this.throttledProviders.set(provider, Date.now() + this.backoffTime);
  }

  _isThrottled(provider) {
    const until = this.throttledProviders.get(provider);
    if (!until) return false;
    if (Date.now() > until) {
      this.throttledProviders.delete(provider);
      return false;
    }
    return true;
  }
}

export const aiManager = new AIManager();
export default aiManager;
