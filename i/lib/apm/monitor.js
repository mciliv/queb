const http = require('http');

const PORT = process.env.APM_PORT || 7243;
const COLORS = { running: '\x1b[32m', idle: '\x1b[90m', error: '\x1b[31m' };

function fetchAgents() {
  return new Promise(resolve => {
    http.get(`http://127.0.0.1:${PORT}/agents`, res => {
      let d = '';
      res.on('data', c => d += c);
      res.on('end', () => { try { resolve(JSON.parse(d)); } catch { resolve([]); } });
    }).on('error', () => resolve([]));
  });
}

async function render() {
  const agents = await fetchAgents();
  process.stdout.write('\x1B[2J\x1B[0f');

  console.log('\x1b[36m' + '='.repeat(60));
  console.log('  AGENT PROCESS MONITOR');
  console.log('='.repeat(60) + '\x1b[0m\n');

  if (!agents.length) {
    console.log('  \x1b[90mNo agents reporting.\x1b[0m');
  } else {
    for (const a of agents) {
      const age = Math.round((Date.now() - a.lastUpdate) / 1000);
      const c = COLORS[a.status] || '\x1b[90m';
      console.log(`${c}[${(a.status || '?').toUpperCase()}]\x1b[0m  ${a.name || a.id}  \x1b[90m${age}s ago\x1b[0m`);
      if (a.agent)   console.log(`  agent: ${a.agent}  provider: ${a.provider || '?'}${a.model ? '/' + a.model : ''}`);
      if (a.prompt)  console.log(`  \x1b[3m"${a.prompt.slice(0, 80)}"\x1b[0m`);
      if (a.error)   console.log(`  \x1b[31m${a.error}\x1b[0m`);
      console.log();
    }
  }

  console.log('\x1b[90m' + '-'.repeat(60));
  console.log(`  ${new Date().toLocaleTimeString()}  |  Ctrl+C to exit`);
  console.log('-'.repeat(60) + '\x1b[0m');
}

if (require.main === module) { setInterval(render, 1000); render(); }
module.exports = { render };
