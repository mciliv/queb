const http = require('http');

const PORT = process.env.APM_PORT || 7243;
const agents = new Map();
const directives = new Map(); // agentId -> [{ message, from, timestamp }]

const server = http.createServer((req, res) => {
  const url = req.url;

  if (req.method === 'POST' && url.startsWith('/ingest/')) {
    const id = url.slice('/ingest/'.length);
    let body = '';
    req.on('data', c => body += c);
    req.on('end', () => {
      try {
        const data = JSON.parse(body);
        agents.set(id, { id, ...agents.get(id), ...data, lastUpdate: Date.now() });
      } catch {}
      res.writeHead(204).end();
    });

  } else if (req.method === 'GET' && url === '/agents') {
    res.writeHead(200, { 'Content-Type': 'application/json' });
    res.end(JSON.stringify(Array.from(agents.values())));

  } else if (req.method === 'POST' && url.startsWith('/directive/')) {
    // Manager sends a directive to a specific agent
    const id = url.slice('/directive/'.length);
    let body = '';
    req.on('data', c => body += c);
    req.on('end', () => {
      try {
        const data = JSON.parse(body);
        if (!directives.has(id)) directives.set(id, []);
        directives.get(id).push({ ...data, timestamp: Date.now() });
      } catch {}
      res.writeHead(204).end();
    });

  } else if (req.method === 'GET' && url.startsWith('/directive/')) {
    // Agent polls for its directives (drains the queue)
    const id = url.slice('/directive/'.length);
    const pending = directives.get(id) || [];
    directives.set(id, []);
    res.writeHead(200, { 'Content-Type': 'application/json' });
    res.end(JSON.stringify(pending));

  } else {
    res.writeHead(404).end();
  }
});

function start() {
  server.listen(PORT, '127.0.0.1', () => {
    console.log(`APM server on http://127.0.0.1:${PORT}`);
  });
}

if (require.main === module) start();
module.exports = { start };
