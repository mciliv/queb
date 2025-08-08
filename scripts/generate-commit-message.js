#!/usr/bin/env node

const { execSync } = require('child_process');
const fs = require('fs');
const path = require('path');

async function generateCommitMessage() {
  // Get git diff
  const diff = execSync('git diff --cached --stat', { encoding: 'utf-8' });
  const detailedDiff = execSync('git diff --cached', { encoding: 'utf-8', maxBuffer: 1024 * 1024 });
  
  // Try OpenAI first
  try {
    const openai = require('openai');
    const apiKey = process.env.OPENAI_API_KEY;
    
    if (apiKey) {
      const client = new openai.OpenAI({ apiKey });
      
      const prompt = `Generate a concise, semantic commit message for these changes:

Summary:
${diff}

First 500 chars of detailed changes:
${detailedDiff.substring(0, 500)}

Guidelines:
- Use conventional commit format (feat:, fix:, refactor:, docs:, style:, test:, chore:)
- Be specific but concise (max 72 chars)
- Focus on WHAT changed and WHY, not HOW
- Single line only`;

      const response = await client.chat.completions.create({
        model: 'gpt-3.5-turbo',
        messages: [{ role: 'user', content: prompt }],
        max_tokens: 100,
        temperature: 0.3
      });
      
      return response.choices[0].message.content.trim();
    }
  } catch (err) {
    console.log('OpenAI unavailable, using fallback...');
  }
  
  // Fallback: analyze diff to generate message
  const files = diff.split('\n').filter(line => line.includes('|'));
  const fileCount = files.length;
  
  if (fileCount === 0) return 'chore: minor updates';
  
  // Analyze file types and changes
  const changes = {
    feat: 0,
    fix: 0,
    test: 0,
    docs: 0,
    style: 0,
    refactor: 0,
    chore: 0
  };
  
  files.forEach(file => {
    const filePath = file.split('|')[0].trim();
    if (filePath.includes('test/') || filePath.includes('.test.')) changes.test++;
    else if (filePath.includes('README') || filePath.includes('.md')) changes.docs++;
    else if (filePath.includes('.css') || filePath.includes('.scss')) changes.style++;
    else if (filePath.includes('package.json') || filePath.includes('config')) changes.chore++;
    else if (detailedDiff.includes('fix') || detailedDiff.includes('bug')) changes.fix++;
    else changes.feat++;
  });
  
  // Determine primary change type
  const primaryType = Object.entries(changes).reduce((a, b) => changes[a[0]] > changes[b[0]] ? a : b)[0];
  
  // Generate message based on analysis
  if (fileCount === 1) {
    const fileName = path.basename(files[0].split('|')[0].trim());
    return `${primaryType}: update ${fileName}`;
  } else {
    const areas = [...new Set(files.map(f => {
      const parts = f.split('|')[0].trim().split('/');
      return parts[0] === 'backend' || parts[0] === 'frontend' ? parts[0] : 'project';
    }))];
    
    return `${primaryType}: update ${areas.join(', ')} (${fileCount} files)`;
  }
}

// Run if called directly
if (require.main === module) {
  generateCommitMessage().then(message => {
    console.log(message);
  }).catch(err => {
    console.error('Error:', err.message);
    console.log('chore: update files');
  });
}

module.exports = { generateCommitMessage };
